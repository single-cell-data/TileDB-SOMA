/**
 * @file   soma_object.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMAObject bindings.
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

void load_soma_object(py::module& m) {
    py::class_<SOMAObject>(m, "SOMAObject")

        .def_static(
            "open",
            [](std::string_view uri,
               OpenMode mode,
               std::shared_ptr<SOMAContext> ctx,
               std::optional<std::pair<uint64_t, uint64_t>> timestamp)
                -> py::object {
                try {
                    auto obj = SOMAObject::open(uri, mode, ctx, timestamp);
                    if (obj->type() == "SOMADataFrame")
                        return py::cast(dynamic_cast<SOMADataFrame&>(*obj));
                    else if (obj->type() == "SOMASparseNDArray")
                        return py::cast(dynamic_cast<SOMASparseNDArray&>(*obj));
                    else if (obj->type() == "SOMADenseNDArray")
                        return py::cast(dynamic_cast<SOMADenseNDArray&>(*obj));
                    else if (obj->type() == "SOMACollection")
                        return py::cast(dynamic_cast<SOMACollection&>(*obj));
                    else if (obj->type() == "SOMAExperiment")
                        return py::cast(dynamic_cast<SOMAExperiment&>(*obj));
                    else if (obj->type() == "SOMAMeasurement")
                        return py::cast(dynamic_cast<SOMAMeasurement&>(*obj));
                    return py::none();
                } catch (...) {
                    return py::none();
                }
            })
        .def_property_readonly("type", &SOMAObject::type);
};
}  // namespace libtiledbsomacpp
