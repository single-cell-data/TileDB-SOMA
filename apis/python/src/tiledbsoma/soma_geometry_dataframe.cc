/**
 * @file   soma_point_cloud_dataframe.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMAGeometryDataFrame bindings.
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

void load_soma_geometry_dataframe(py::module& m) {
    py::class_<SOMAGeometryDataFrame, SOMAArray, SOMAObject>(
        m, "SOMAGeometryDataFrame")

        .def_static(
            "create",
            [](std::string_view uri,
               py::object py_schema,
               py::object index_column_info,
               std::vector<std::string> axis_names,
               std::vector<std::optional<std::string>> axis_units,
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

                SOMACoordinateSpace coord_space{axis_names, axis_units};

                try {
                    SOMAGeometryDataFrame::create(
                        uri,
                        std::make_unique<ArrowSchema>(schema),
                        ArrowTable(
                            std::make_unique<ArrowArray>(index_column_array),
                            std::make_unique<ArrowSchema>(index_column_schema)),
                        coord_space,
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
            "axis_names"_a,
            "axis_units"_a,
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
                &SOMAGeometryDataFrame::open),
            "uri"_a,
            "mode"_a,
            "context"_a,
            py::kw_only(),
            "timestamp"_a = py::none(),
            py::call_guard<py::gil_scoped_release>())

        .def_static("exists", &SOMAGeometryDataFrame::exists)
        .def_property_readonly(
            "index_column_names", &SOMAGeometryDataFrame::index_column_names)
        .def_property_readonly(
            "count",
            &SOMAGeometryDataFrame::count,
            py::call_guard<py::gil_scoped_release>());
}
}  // namespace libtiledbsomacpp
