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

void load_column_filter(py::module& m) {
    py::class_<CoordinateValueFilter>(m, "CoordinateValueFilter")
        .def(py::init<const SOMAArray&>())
        .def(
            "add_arrow_points", []() {}, "column_index"_a, "points"_a)
        .def("add_int64_slice", &CoordinateValueFilter::<int64_t>add_slice)
        .def("add_int64_point_selection", &SOMAIndexColumnFilter::<int64_t>add_points)
        .def("add_point_selection", &CoordinateValueFilter::add_points)
        .def("add_slice", &CoordinateValueFilter::add_slice)
        .def("add_string_slice", &CoordinateValueFilter::add_slice)
        .def("add_string_point_selection", &CoordinateValueFilter::add_string_points);
}

}  // namespace libtiledbsomacpp
