/**
 * @file   soma_vfs.cc
 *
 * @section LICENSE
 *
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
                    TPY_ERROR_LOC(
                        std::format("URI {} is not a valid URI", uri));
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
