/**
 * @file   soma_collection.cc
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
 * This file defines the SOMACollection bindings.
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

void load_soma_collection(py::module& m) {
    py::class_<SOMACollection, SOMAGroup, SOMAObject>(m, "SOMACollection")
        .def_static(
            "open",
            py::overload_cast<
                std::string_view,
                OpenMode,
                std::shared_ptr<SOMAContext>,
                std::optional<std::pair<uint64_t, uint64_t>>>(
                &SOMACollection::open),
            "uri"_a,
            py::kw_only(),
            "mode"_a,
            "context"_a,
            "timestamp"_a = py::none())
        .def(
            "__iter__",
            [](SOMACollection& collection) {
                return py::make_iterator(collection.begin(), collection.end());
            },
            py::keep_alive<0, 1>())
        .def("get", &SOMACollection::get);

    py::class_<SOMAExperiment, SOMACollection, SOMAGroup, SOMAObject>(
        m, "SOMAExperiment");

    py::class_<SOMAMeasurement, SOMACollection, SOMAGroup, SOMAObject>(
        m, "SOMAMeasurement");

    py::class_<SOMAScene, SOMACollection, SOMAGroup, SOMAObject>(
        m, "SOMAScene");

    py::class_<SOMAMultiscaleImage, SOMACollection, SOMAGroup, SOMAObject>(
        m, "SOMAMultiscaleImage");
}
}  // namespace libtiledbsomacpp
