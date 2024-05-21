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
               std::shared_ptr<SOMAContext> context,
               std::optional<std::pair<uint64_t, uint64_t>> timestamp,
               std::optional<std::string> clib_type) -> py::object {
                try {
                    auto soma_obj = SOMAObject::open(
                        uri, mode, context, timestamp, clib_type);
                    auto soma_obj_type = soma_obj->type();

                    if (soma_obj_type.has_value()) {
                        std::transform(
                            soma_obj_type->begin(),
                            soma_obj_type->end(),
                            soma_obj_type->begin(),
                            [](unsigned char c) { return std::tolower(c); });
                    }

                    if (soma_obj_type == "somadataframe")
                        return py::cast(
                            dynamic_cast<SOMADataFrame&>(*soma_obj));
                    else if (soma_obj_type == "somasparsendarray")
                        return py::cast(
                            dynamic_cast<SOMASparseNDArray&>(*soma_obj));
                    else if (soma_obj_type == "somadensendarray")
                        return py::cast(
                            dynamic_cast<SOMADenseNDArray&>(*soma_obj));
                    else if (soma_obj_type == "somacollection")
                        return py::cast(
                            dynamic_cast<SOMACollection&>(*soma_obj));
                    else if (soma_obj_type == "somaexperiment")
                        return py::cast(
                            dynamic_cast<SOMAExperiment&>(*soma_obj));
                    else if (soma_obj_type == "somameasurement")
                        return py::cast(
                            dynamic_cast<SOMAMeasurement&>(*soma_obj));
                    return py::none();
                } catch (...) {
                    return py::none();
                }
            },
            "uri"_a,
            "mode"_a,
            "context"_a,
            py::kw_only(),
            "timestamp"_a = py::none(),
            "clib_type"_a = py::none())
        .def_property_readonly("type", &SOMAObject::type);
};
}  // namespace libtiledbsomacpp
