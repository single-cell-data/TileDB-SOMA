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
        .def_static(
            "create",
            [](std::shared_ptr<SOMAContext> ctx,
               std::string_view uri,
               std::string soma_type,
               std::optional<TimestampRange> timestamp) {
                try {
                    SOMAGroup::create(ctx, uri, soma_type, timestamp);
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            },
            py::kw_only(),
            "ctx"_a,
            "uri"_a,
            "soma_type"_a,
            "timestamp"_a = py::none())
        .def("__enter__", [](SOMAGroup& group) { return group; })
        .def(
            "__exit__",
            [](SOMAGroup& group,
               py::object exc_type,
               py::object exc_value,
               py::object traceback) { group.close(); })
        .def_property_readonly(
            "mode",
            [](SOMAGroup& group) {
                return group.mode() == OpenMode::read ? "r" : "w";
            })
        .def("reopen", &SOMAGroup::reopen)
        .def("close", &SOMAGroup::close)
        .def_property_readonly(
            "closed",
            [](SOMAGroup& group) -> bool { return not group.is_open(); })
        .def_property_readonly("uri", &SOMAGroup::uri)
        .def("context", &SOMAGroup::ctx)
        .def("has", &SOMAGroup::has)
        .def(
            "add",
            &SOMAGroup::set,
            "uri"_a,
            "uri_type"_a,
            "name"_a,
            "soma_type"_a)
        .def("count", &SOMAGroup::count)
        .def("remove", &SOMAGroup::del)
        .def("members", &SOMAGroup::members_map)
        .def_property_readonly(
            "timestamp",
            [](SOMAGroup& group) -> py::object {
                if (!group.timestamp().has_value())
                    return py::none();
                return py::cast(group.timestamp()->second);
            })
        .def_property_readonly(
            "meta",
            [](SOMAGroup& group) -> py::dict {
                return meta(group.get_metadata());
            })
        .def(
            "set_metadata",
            [](SOMAGroup& group, const std::string& key, py::array value) {
                set_metadata(group, key, value);
            })
        .def("delete_metadata", &SOMAGroup::delete_metadata)
        .def("has_metadata", &SOMAGroup::has_metadata)
        .def("metadata_num", &SOMAGroup::metadata_num);
}
}  // namespace libtiledbsomacpp
