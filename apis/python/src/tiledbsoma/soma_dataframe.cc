/**
 * @file   soma_dataframe.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMADataFrame bindings.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_soma_dataframe(py::module& m) {
    py::class_<SOMADataFrame, SOMAArray, SOMAObject>(m, "SOMADataFrame")

        .def_static(
            "create",
            [](std::string_view uri,
               py::object py_schema,
               py::object index_column_info,
               std::shared_ptr<SOMAContext> context,
               PlatformConfig platform_config,
               std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                ArrowSchema schema;
                uintptr_t schema_ptr = (uintptr_t)(&schema);
                py_schema.attr("_export_to_c")(schema_ptr);

                // Please see
                // https://github.com/single-cell-data/TileDB-SOMA/issues/2869
                // for the reasoning here.
                //
                // TL;DR:
                // * The table has an `ArrowSchema`; each of its children
                //   is also an `ArrowSchema`.
                // * Arrow fields are nullable by default in the user API.
                // * There is a field-level nullability flag, _and_ users
                //   can set a "nullable" metadata as well.
                // * In the absence of metadata, respect the flag we get.
                // * In the present of metdata with "nullable", let that
                //   override.

                auto metadata = py_schema.attr("metadata");
                if (py::hasattr(metadata, "get")) {
                    for (int64_t i = 0; i < schema.n_children; ++i) {
                        auto child = schema.children[i];
                        auto val = metadata.attr("get")(
                            py::str(child->name).attr("encode")("utf-8"));

                        if (!val.is(py::none()) &&
                            val.cast<std::string>() == "nullable") {
                            child->flags |= ARROW_FLAG_NULLABLE;
                        }
                    }
                }

                ArrowSchema index_column_schema;
                ArrowArray index_column_array;
                uintptr_t
                    index_column_schema_ptr = (uintptr_t)(&index_column_schema);
                uintptr_t
                    index_column_array_ptr = (uintptr_t)(&index_column_array);
                index_column_info.attr("_export_to_c")(
                    index_column_array_ptr, index_column_schema_ptr);

                try {
                    SOMADataFrame::create(
                        uri,
                        std::make_unique<ArrowSchema>(schema),
                        ArrowTable(
                            std::make_unique<ArrowArray>(index_column_array),
                            std::make_unique<ArrowSchema>(index_column_schema)),
                        context,
                        platform_config,
                        timestamp);
                } catch (const std::out_of_range& e) {
                    throw py::type_error(e.what());
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
                schema.release(&schema);
            },
            "uri"_a,
            py::kw_only(),
            "schema"_a,
            "index_column_info"_a,
            "ctx"_a,
            "platform_config"_a,
            "timestamp"_a = py::none())

        .def_static(
            "open",
            py::overload_cast<
                std::string_view,
                OpenMode,
                std::shared_ptr<SOMAContext>,
                std::optional<std::pair<uint64_t, uint64_t>>>(
                &SOMADataFrame::open),
            "uri"_a,
            "mode"_a,
            "context"_a,
            py::kw_only(),
            "timestamp"_a = py::none(),
            py::call_guard<py::gil_scoped_release>())

        .def_static("exists", &SOMADataFrame::exists)
        .def_property_readonly(
            "index_column_names", &SOMADataFrame::index_column_names)

        .def_property_readonly(
            "count",
            &SOMADataFrame::count,
            py::call_guard<py::gil_scoped_release>())
        .def_property_readonly(
            "maybe_soma_joinid_shape", &SOMADataFrame::maybe_soma_joinid_shape)
        .def_property_readonly(
            "maybe_soma_joinid_maxshape",
            &SOMADataFrame::maybe_soma_joinid_maxshape)
        .def_property_readonly(
            "tiledbsoma_has_upgraded_domain", &SOMAArray::has_current_domain)

        .def(
            "resize_soma_joinid_shape",
            [](SOMADataFrame& sdf,
               int64_t newshape,
               std::string function_name_for_messages) {
                try {
                    sdf.resize_soma_joinid_shape(
                        newshape, function_name_for_messages);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a,
            "function_name_for_messages"_a)

        .def(
            "can_resize_soma_joinid_shape",
            [](SOMADataFrame& sdf,
               int64_t newshape,
               std::string function_name_for_messages) {
                try {
                    return sdf.can_resize_soma_joinid_shape(
                        newshape, function_name_for_messages);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a,
            "function_name_for_messages"_a)

        .def(
            "upgrade_soma_joinid_shape",
            [](SOMADataFrame& sdf,
               int64_t newshape,
               std::string function_name_for_messages) {
                try {
                    sdf.upgrade_soma_joinid_shape(
                        newshape, function_name_for_messages);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a,
            "function_name_for_messages"_a)

        .def(
            "can_upgrade_soma_joinid_shape",
            [](SOMADataFrame& sdf,
               int64_t newshape,
               std::string function_name_for_messages) {
                try {
                    return sdf.can_upgrade_soma_joinid_shape(
                        newshape, function_name_for_messages);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a,
            "function_name_for_messages"_a)

        .def(
            "upgrade_domain",
            [](SOMADataFrame& sdf,
               py::object pyarrow_domain_table,
               std::string function_name_for_messages) {
                ArrowArray pyarrow_domain_array;
                ArrowSchema pyarrow_domain_schema;
                uintptr_t nanoarrow_domain_array_ptr =
                    (uintptr_t)(&pyarrow_domain_array);
                uintptr_t nanoarrow_domain_schema_ptr =
                    (uintptr_t)(&pyarrow_domain_schema);
                pyarrow_domain_table.attr("_export_to_c")(
                    nanoarrow_domain_array_ptr, nanoarrow_domain_schema_ptr);
                ArrowTable nanoarrow_domain_table(
                    std::make_unique<ArrowArray>(pyarrow_domain_array),
                    std::make_unique<ArrowSchema>(pyarrow_domain_schema));
                try {
                    sdf.upgrade_domain(
                        nanoarrow_domain_table, function_name_for_messages);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "pyarrow_domain_table"_a,
            "function_name_for_messages"_a)

        .def(
            "can_upgrade_domain",
            [](SOMADataFrame& sdf,
               py::object pyarrow_domain_table,
               std::string function_name_for_messages) {
                ArrowArray pyarrow_domain_array;
                ArrowSchema pyarrow_domain_schema;
                uintptr_t nanoarrow_domain_array_ptr =
                    (uintptr_t)(&pyarrow_domain_array);
                uintptr_t nanoarrow_domain_schema_ptr =
                    (uintptr_t)(&pyarrow_domain_schema);
                pyarrow_domain_table.attr("_export_to_c")(
                    nanoarrow_domain_array_ptr, nanoarrow_domain_schema_ptr);
                ArrowTable nanoarrow_domain_table(
                    std::make_unique<ArrowArray>(pyarrow_domain_array),
                    std::make_unique<ArrowSchema>(pyarrow_domain_schema));
                try {
                    return sdf.can_upgrade_domain(
                        nanoarrow_domain_table, function_name_for_messages);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "pyarrow_domain_table"_a,
            "function_name_for_messages"_a)

        .def(
            "change_domain",
            [](SOMADataFrame& sdf,
               py::object pyarrow_domain_table,
               std::string function_name_for_messages) {
                ArrowArray pyarrow_domain_array;
                ArrowSchema pyarrow_domain_schema;
                uintptr_t nanoarrow_domain_array_ptr =
                    (uintptr_t)(&pyarrow_domain_array);
                uintptr_t nanoarrow_domain_schema_ptr =
                    (uintptr_t)(&pyarrow_domain_schema);
                pyarrow_domain_table.attr("_export_to_c")(
                    nanoarrow_domain_array_ptr, nanoarrow_domain_schema_ptr);
                ArrowTable nanoarrow_domain_table(
                    std::make_unique<ArrowArray>(pyarrow_domain_array),
                    std::make_unique<ArrowSchema>(pyarrow_domain_schema));
                try {
                    sdf.change_domain(
                        nanoarrow_domain_table, function_name_for_messages);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "pyarrow_domain_table"_a,
            "function_name_for_messages"_a)

        .def(
            "can_change_domain",
            [](SOMADataFrame& sdf,
               py::object pyarrow_domain_table,
               std::string function_name_for_messages) {
                ArrowArray pyarrow_domain_array;
                ArrowSchema pyarrow_domain_schema;
                uintptr_t nanoarrow_domain_array_ptr =
                    (uintptr_t)(&pyarrow_domain_array);
                uintptr_t nanoarrow_domain_schema_ptr =
                    (uintptr_t)(&pyarrow_domain_schema);
                pyarrow_domain_table.attr("_export_to_c")(
                    nanoarrow_domain_array_ptr, nanoarrow_domain_schema_ptr);
                ArrowTable nanoarrow_domain_table(
                    std::make_unique<ArrowArray>(pyarrow_domain_array),
                    std::make_unique<ArrowSchema>(pyarrow_domain_schema));
                try {
                    return sdf.can_change_domain(
                        nanoarrow_domain_table, function_name_for_messages);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "pyarrow_domain_table"_a,
            "function_name_for_messages"_a)

        .def(
            "_update_dataframe_schema",
            &SOMADataFrame::update_dataframe_schema);
}

}  // namespace libtiledbsomacpp
