/**
 * @file   soma_context.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMAContext bindings.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <memory>
#include <tiledbsoma/tiledbsoma>

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_soma_context(py::module& m) {
    py::class_<SOMAContext, std::shared_ptr<SOMAContext>>(m, "SOMAContext")
        .def(py::init<>())
        .def(py::init<std::map<std::string, std::string>>())
        .def("config", &SOMAContext::tiledb_config)
        .def("tiledb_ctx", &SOMAContext::tiledb_ctx);
};
}  // namespace libtiledbsomacpp
