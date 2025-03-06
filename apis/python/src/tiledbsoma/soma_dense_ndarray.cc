/**
 * @file   soma_dense_ndarray.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMADenseNDArray bindings.
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

void write(SOMAArray& array, py::array data) {
    ManagedQuery mq(array, array.ctx()->tiledb_ctx());

    py::buffer_info data_info = data.request();
    mq.setup_write_column(
        "soma_data",
        data.size(),
        (const void*)data_info.ptr,
        (uint64_t*)nullptr);

    try {
        mq.submit_write();
        mq.close();
    } catch (const std::exception& e) {
        TPY_ERROR_LOC(e.what());
    }
}

void load_soma_dense_ndarray(py::module& m) {
    py::class_<SOMADenseNDArray, SOMAArray, SOMAObject>(m, "SOMADenseNDArray")

        .def_static(
            "create",
            [](std::string_view uri,
               std::string format,
               py::object index_column_info,
               std::shared_ptr<SOMAContext> context,
               PlatformConfig platform_config,
               std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                ArrowSchema index_column_schema;
                ArrowArray index_column_array;
                uintptr_t
                    index_column_schema_ptr = (uintptr_t)(&index_column_schema);
                uintptr_t
                    index_column_array_ptr = (uintptr_t)(&index_column_array);
                index_column_info.attr("_export_to_c")(
                    index_column_array_ptr, index_column_schema_ptr);

                try {
                    SOMADenseNDArray::create(
                        uri,
                        format,
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

                index_column_array.release(&index_column_array);
                index_column_schema.release(&index_column_schema);
            },
            "uri"_a,
            py::kw_only(),
            "format"_a,
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
                &SOMADenseNDArray::open),
            "uri"_a,
            "mode"_a,
            "context"_a,
            py::kw_only(),
            "timestamp"_a = py::none(),
            py::call_guard<py::gil_scoped_release>())

        .def_static("exists", &SOMADenseNDArray::exists)

        .def("write", write)

        .def_property_readonly("shape", &SOMADenseNDArray::shape)
        .def_property_readonly("maxshape", &SOMADenseNDArray::maxshape);
}
}  // namespace libtiledbsomacpp
