/**
 * @file   soma_vfs.cc
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
 * This file defines the VFS bindings.
 */

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;
using VFSFilebuf = tiledb::impl::VFSFilebuf;

// TODO This temporary workaround prevents namespace clash with tiledb-py.
// Bind tiledb::VFS directly once tiledb-py dependency is removed
class SOMAVFS : public tiledb::VFS {
   public:
    using tiledb::VFS::VFS;
};

void load_soma_vfs(py::module& m) {
    py::class_<SOMAVFS>(m, "SOMAVFS")
        .def(
            py::init([](std::shared_ptr<SOMAContext> context) {
                return SOMAVFS(*context->tiledb_ctx());
            }),
            "ctx"_a);

    py::class_<VFSFilebuf>(m, "SOMAVFSFilebuf")
        .def(py::init<const SOMAVFS&>())
        .def(
            "open",
            [](VFSFilebuf& buf, const std::string& uri) {
                return buf.open(uri, std::ios::in);
            })
        .def("close", &VFSFilebuf::close, "should_throw"_a = true);
}
}  // namespace libtiledbsomacpp
