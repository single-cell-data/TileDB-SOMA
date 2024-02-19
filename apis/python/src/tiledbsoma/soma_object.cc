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

void load_soma_object(py::module &m) {
    py::class_<SOMAObject>(m, "SOMAObject")

    .def_static("open", [](std::string uri, 
                           OpenMode mode, 
                           std::map<std::string, std::string> config, 
                           std::optional<std::pair<uint64_t, uint64_t>> timestamp) -> py::object {
        if(mode == OpenMode::write)
            TPY_ERROR_LOC("SOMAObjects for write mode not handled in Python API yet.");

        try{
            auto obj = SOMAObject::open(uri, mode, config, timestamp);
            if (obj->type() == "SOMADataFrame")
                return py::cast(dynamic_cast<SOMADataFrame&>(*obj));
        }
        catch(...){
            TPY_ERROR_LOC("SOMAObject not handled in Python API yet.");
        }
        })

        .def_property_readonly("type", &SOMAObject::type);
    };
}

