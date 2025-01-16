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

// TODO This temporary workaround prevents namespace clash with tiledb-py.
// Bind tiledb::VFS directly once tiledb-py dependency is removed
class SOMAVFS : public tiledb::VFS {
   public:
    using tiledb::VFS::VFS;
    SOMAVFS(const tiledb::VFS& vfs)
        : tiledb::VFS(vfs) {
    }
};

class SOMAVFSFilebuf : public tiledb::impl::VFSFilebuf {
   private:
    std::streamsize offset_ = 0;
    SOMAVFS vfs_;

   public:
    SOMAVFSFilebuf(const VFS& vfs)
        : tiledb::impl::VFSFilebuf(vfs)
        , vfs_(vfs){};

    std::streamsize seek(std::streamsize offset, uint64_t whence) {
        if (whence == 0) {
            offset_ = seekoff(offset, std::ios::beg, std::ios::in);
        } else if (whence == 1) {
            offset_ += seekoff(offset, std::ios::cur, std::ios::in);
        } else if (whence == 2) {
            offset_ = vfs_.file_size(get_uri()) -
                      seekoff(offset, std::ios::end, std::ios::in);
        } else {
            TPY_ERROR_LOC(
                "whence must be equal to SEEK_SET, SEEK_CUR, SEEK_END");
        }

        return offset_;
    }

    py::bytes read(std::streamsize size) {
        int64_t nbytes = (size < 0 || size > showmanyc()) ? showmanyc() : size;
        if (nbytes <= 0) {
            return py::bytes("");
        }

        py::gil_scoped_release release;
        std::string buffer(nbytes, '\0');
        offset_ += xsgetn(&buffer[0], nbytes);
        py::gil_scoped_acquire acquire;

        return py::bytes(buffer);
    }

    std::streamsize tell() {
        return offset_;
    }
};

void load_soma_vfs(py::module& m) {
    py::class_<SOMAVFS>(m, "SOMAVFS")
        .def(
            py::init([](std::shared_ptr<SOMAContext> context) {
                return SOMAVFS(*context->tiledb_ctx());
            }),
            "ctx"_a);

    py::class_<SOMAVFSFilebuf>(m, "SOMAVFSFilebuf")
        .def(py::init<const SOMAVFS&>())
        .def(
            "open",
            [](SOMAVFSFilebuf& buf, const std::string& uri) {
                auto fb = buf.open(uri, std::ios::in);
                if (fb == nullptr) {
                    // No std::format in C++17, and fmt::format is overkill here
                    std::stringstream ss;
                    ss << "URI " << uri << " is not a valid URI";
                    TPY_ERROR_LOC(ss.str());
                }
                return fb;
            },
            py::call_guard<py::gil_scoped_release>())
        .def("read", &SOMAVFSFilebuf::read, "size"_a = -1)
        .def("tell", &SOMAVFSFilebuf::tell)
        .def(
            "seek",
            &SOMAVFSFilebuf::seek,
            "offset"_a,
            "whence"_a = 0,
            py::call_guard<py::gil_scoped_release>())
        .def("close", &SOMAVFSFilebuf::close, "should_throw"_a = true);
}
}  // namespace libtiledbsomacpp
