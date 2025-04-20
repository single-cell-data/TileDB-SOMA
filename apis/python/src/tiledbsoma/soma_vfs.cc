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

#include <streambuf>

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

class SOMAFileHandle {
   public:
    class SOMAFilebuf : public tiledb::impl::VFSFilebuf {
       public:
        using tiledb::impl::VFSFilebuf::seekoff;
        using tiledb::impl::VFSFilebuf::showmanyc;
        using tiledb::impl::VFSFilebuf::VFSFilebuf;
        using tiledb::impl::VFSFilebuf::xsgetn;
    };

   private:
    std::streamsize offset_ = 0;
    std::ios::openmode openmode_ = std::ios::in;  // hardwired to read-only

    // Ensure proper order of destruction
    SOMAFilebuf buf_;
    tiledb::VFS vfs_;
    std::shared_ptr<SOMAContext> ctx_;

   public:
    SOMAFileHandle(const std::string& uri, std::shared_ptr<SOMAContext> ctx)
        : vfs_(tiledb::VFS(*ctx->tiledb_ctx()))
        , buf_(SOMAFilebuf(vfs_))
        , ctx_(ctx) {
        if (buf_.open(uri, openmode_) == nullptr) {
            // No std::format in C++17, and fmt::format is overkill here
            std::stringstream ss;
            ss << "URI " << uri << " is not a valid URI";
            TPY_ERROR_LOC(ss.str());
        }
    }

    std::streamsize seek(std::streamsize offset, uint64_t whence) {
        if (closed()) {
            TPY_ERROR_LOC("File must be open before performing seek");
        }

        if (whence == 0) {
            offset_ = buf_.seekoff(offset, std::ios::beg, std::ios::in);
        } else if (whence == 1) {
            offset_ += buf_.seekoff(offset, std::ios::cur, std::ios::in);
        } else if (whence == 2) {
            offset_ = vfs_.file_size(buf_.get_uri()) -
                      buf_.seekoff(offset, std::ios::end, std::ios::in);
        } else {
            TPY_ERROR_LOC(
                "whence must be equal to SEEK_SET, SEEK_CUR, SEEK_END");
        }

        return offset_;
    }

    py::bytes read(std::streamsize size) {
        if (closed()) {
            TPY_ERROR_LOC("File must be open before performing read");
        }

        int64_t nbytes = (size < 0 || size > buf_.showmanyc()) ?
                             buf_.showmanyc() :
                             size;
        if (nbytes <= 0) {
            return py::bytes("");
        }

        py::gil_scoped_release release;
        std::string buffer(nbytes, '\0');
        offset_ += buf_.xsgetn(&buffer[0], nbytes);
        py::gil_scoped_acquire acquire;

        return py::bytes(buffer);
    }

    std::streamsize readinto(py::buffer buffer) {
        if (closed()) {
            TPY_ERROR_LOC("File must be open before performing readinto");
        }

        py::buffer_info info = buffer.request();
        if (info.ndim != 1)
            throw std::runtime_error("Expected a 1-dimensional byte array");
        if (info.readonly)
            throw std::runtime_error("Cannot write to a read-only buffer");

        auto nbytes = info.size;
        if (nbytes <= 0)
            return 0;

        py::gil_scoped_release release;
        auto bytes_read = buf_.xsgetn(static_cast<char*>(info.ptr), nbytes);
        offset_ += bytes_read;
        py::gil_scoped_acquire acquire;

        return bytes_read;
    }

    std::streamsize tell() {
        return offset_;
    }

    void close(bool should_throw) {
        buf_.close(should_throw);
    }

    bool readable() {
        return (openmode_ & std::ios::in) != 0;
    }

    bool writable() {
        return (openmode_ & std::ios::out) != 0;
    }

    bool closed() {
        return !buf_.is_open();
    }

    bool seekable() {
        return true;
    }
};

void load_soma_vfs(py::module& m) {
    py::class_<SOMAFileHandle>(m, "SOMAFileHandle")
        .def(
            py::init<const std::string&, std::shared_ptr<SOMAContext>>(),
            "uri"_a,
            "ctx"_a,
            py::call_guard<py::gil_scoped_release>())
        .def("read", &SOMAFileHandle::read, "size"_a = -1)
        .def("readinto", &SOMAFileHandle::readinto, "buffer"_a)
        .def("flush", [](SOMAFileHandle& buf) {})
        .def("tell", &SOMAFileHandle::tell)
        .def("readable", &SOMAFileHandle::readable)
        .def("writable", &SOMAFileHandle::writable)
        .def_property_readonly("closed", &SOMAFileHandle::closed)
        .def("seekable", &SOMAFileHandle::seekable)
        .def(
            "seek",
            &SOMAFileHandle::seek,
            "offset"_a,
            "whence"_a = 0,
            py::call_guard<py::gil_scoped_release>())
        .def("close", &SOMAFileHandle::close, "should_throw"_a = true);
}
}  // namespace libtiledbsomacpp
