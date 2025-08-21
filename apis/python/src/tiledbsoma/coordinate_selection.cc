/**
 * @file   soma_column_filter.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines column filters for querying SOMA.
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

template <typename T>
void add_points(CoordinateValueFilters& value_filter, int64_t col_index, const std::vector<T>& values) {
    value_filter.add_points<T>(col_index, PointSelection<T>(values));
}

template <typename T>
void add_slice(CoordinateValueFilters& value_filter, int64_t col_index, std::optional<T> start, std::optional<T> stop) {
    value_filter.add_slice<T>(col_index, SliceSelection<T>(start, stop));
}

void load_coordinate_selection(py::module& m) {
    py::class_<CoordinateValueFilters>(m, "CoordinateValueFilters")
        .def(py::init([](SOMAArray array) { return array.create_coordinate_value_filter(); }), py::arg("array"))
        .def(
            "add_arrow_points",
            [](CoordinateValueFilters& value_filter, int64_t col_index, py::object py_arrow_array) {
                // Create a list of array chunks
                py::list array_chunks;
                if (py::hasattr(py_arrow_array, "chunks")) {
                    array_chunks = py_arrow_array.attr("chunks").cast<py::list>();
                } else {
                    array_chunks.append(py_arrow_array);
                }

                for (const pybind11::handle array_handle : array_chunks) {
                    ArrowSchema arrow_schema;
                    ArrowArray arrow_array;
                    uintptr_t arrow_schema_ptr = (uintptr_t)(&arrow_schema);
                    uintptr_t arrow_array_ptr = (uintptr_t)(&arrow_array);

                    // Call handle._export_to_c to get arrow array and schema
                    //
                    // If ever a NumPy array gets in here, there will be an
                    // exception like "AttributeError: 'numpy.ndarray' object
                    // has no attribute '_export_to_c'".
                    array_handle.attr("_export_to_c")(arrow_array_ptr, arrow_schema_ptr);

                    auto coords = array_handle.attr("tolist")();

                    try {
                        if (!strcmp(arrow_schema.format, "l")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<int64_t>>());
                        } else if (!strcmp(arrow_schema.format, "i")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<int32_t>>());
                        } else if (!strcmp(arrow_schema.format, "s")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<int16_t>>());
                        } else if (!strcmp(arrow_schema.format, "c")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<int8_t>>());
                        } else if (!strcmp(arrow_schema.format, "L")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<uint64_t>>());
                        } else if (!strcmp(arrow_schema.format, "I")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<uint32_t>>());
                        } else if (!strcmp(arrow_schema.format, "S")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<uint16_t>>());
                        } else if (!strcmp(arrow_schema.format, "C")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<uint8_t>>());
                        } else if (!strcmp(arrow_schema.format, "f")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<float_t>>());
                        } else if (!strcmp(arrow_schema.format, "g")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<double_t>>());
                        } else if (
                            !strcmp(arrow_schema.format, "u") || !strcmp(arrow_schema.format, "U") ||
                            !strcmp(arrow_schema.format, "z") || !strcmp(arrow_schema.format, "Z")) {
                            add_points(value_filter, col_index, coords.cast<std::vector<std::string>>());
                        } else if (
                            !strcmp(arrow_schema.format, "tss:") || !strcmp(arrow_schema.format, "tsm:") ||
                            !strcmp(arrow_schema.format, "tsu:") || !strcmp(arrow_schema.format, "tsn:")) {
                            // convert the Arrow Array to int64
                            auto pa = py::module::import("pyarrow");
                            coords = array_handle.attr("cast")(pa.attr("int64")()).attr("tolist")();
                            add_points(value_filter, col_index, coords.cast<std::vector<int64_t>>());
                        } else {
                            TPY_ERROR_LOC(
                                "[pytiledbsoma] set_dim_points: type={} not "
                                "supported" +
                                std::string(arrow_schema.format));
                        }
                    } catch (const std::exception& e) {
                        throw TileDBSOMAError(e.what());
                    }

                    // Release arrow schema
                    arrow_schema.release(&arrow_schema);
                }
            })
        .def("add_slice_int8", add_slice<int8_t>)
        .def("add_slice_int16", add_slice<int16_t>)
        .def("add_slice_int32", add_slice<int32_t>)
        .def("add_slice_int64", add_slice<int64_t>)
        .def("add_slice_uint8", add_slice<uint8_t>)
        .def("add_slice_uint16", add_slice<uint16_t>)
        .def("add_slice_uint32", add_slice<uint32_t>)
        .def("add_slice_uint64", add_slice<uint64_t>)
        .def("add_slice_float", add_slice<float_t>)
        .def("add_slice_double", add_slice<double_t>)
        .def("add_slice_string", add_slice<std::string>)
        .def("add_points_int8", add_points<int8_t>)
        .def("add_points_int16", add_points<int16_t>)
        .def("add_points_int32", add_points<int32_t>)
        .def("add_points_int64", add_points<int64_t>)
        .def("add_points_uint8", add_points<uint8_t>)
        .def("add_points_uint16", add_points<uint16_t>)
        .def("add_points_uint32", add_points<uint32_t>)
        .def("add_points_uint64", add_points<uint64_t>)
        .def("add_points_float", add_points<float_t>)
        .def("add_points_double", add_points<double_t>)
        .def("add_points_string", add_points<std::string>);
}

}  // namespace libtiledbsomacpp
