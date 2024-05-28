/**
 * @file   soma_group.cc
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
 * THE SOFTWARE.Os
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMAGroup bindings.
 */

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_soma_group(py::module& m) {
    py::class_<SOMAGroup, SOMAObject>(m, "SOMAGroup")
        .def_property_readonly(
            "mode",
            [](SOMAGroup& group) {
                return group.mode() == OpenMode::read ? "r" : "w";
            })
        .def("close", &SOMAGroup::close)
        .def_property_readonly(
            "closed",
            [](SOMAGroup& group) -> bool { return not group.is_open(); })
        .def_property_readonly("uri", &SOMAGroup::uri)
        .def("context", &SOMAGroup::ctx)
        .def("has", &SOMAGroup::has)
        .def("count", &SOMAGroup::count)
        .def("member_to_uri_mapping", &SOMAGroup::member_to_uri_mapping)
        .def("timestamp", &SOMAGroup::timestamp)
        .def_property_readonly(
            "meta",
            [](SOMAGroup& group) -> py::dict {
                return meta(group.get_metadata());
            })
        .def("set_metadata", &SOMAGroup::set_metadata)
        .def("delete_metadata", &SOMAGroup::delete_metadata)
        .def("has_metadata", &SOMAGroup::has_metadata)
        .def("metadata_num", &SOMAGroup::metadata_num);
}
}  // namespace libtiledbsomacpp
