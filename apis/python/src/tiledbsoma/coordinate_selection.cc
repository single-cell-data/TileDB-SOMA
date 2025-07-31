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
void add_points(CoordinateValueFilter& value_filter, int64_t col_index, const std::vector<T>& values) {
    value_filter.add_points<T>(col_index, SOMAPointSelection<T>(values));
}

template <typename T>
void add_slice(CoordinateValueFilter& value_filter, int64_t col_index, std::optional<T> start, std::optional<T> stop) {
    value_filter.add_slice<T>(col_index, SOMASliceSelection<T>(start, stop));
}

void load_coordinate_selection(py::module& m) {
    py::class_<CoordinateValueFilter>(m, "CoordinateValueFilter")
        .def(py::init([](SOMAArray array) { return array.create_coordinate_value_filter(); }), py::arg("array"))
        .def(
            "add_slice_string",
            [](CoordinateValueFilter& value_filter,
               int64_t col_index,
               const std::optional<std::string>& start,
               const std::optional<std::string>& stop) {
                value_filter.add_slice(col_index, SOMASliceSelection<std::string>(start, stop));
            })
        .def("add_slice_int8", add_slice<int8_t>)
        .def("add_slice_int16", add_slice<int16_t>)
        .def("add_slice_int32", add_slice<int32_t>)
        .def("add_slice_int64", add_slice<int64_t>)
        .def("add_slice_uint8", add_slice<uint8_t>)
        .def("add_slice_uint16", add_slice<uint16_t>)
        .def("add_slice_uint32", add_slice<uint32_t>)
        .def("add_slice_uint64", add_slice<uint64_t>)
        .def(
            "add_point_string",
            [](CoordinateValueFilter& value_filter, int64_t col_index, const std::vector<std::string>& values) {
                value_filter.add_points(col_index, SOMAPointSelection<std::string>(values));
            })
        .def("add_point_int8", add_points<int8_t>)
        .def("add_point_int16", add_points<int16_t>)
        .def("add_point_int32", add_points<int32_t>)
        .def("add_point_int64", add_points<int64_t>)
        .def("add_point_uint8", add_points<uint8_t>)
        .def("add_point_uint16", add_points<uint16_t>)
        .def("add_point_uint32", add_points<uint32_t>)
        .def("add_point_uint64", add_points<uint64_t>);
}

}  // namespace libtiledbsomacpp
