#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

using namespace tiledbsoma;

namespace py = pybind11;

namespace tiledbsoma {
void init_soma_dense_ndarray(py::module &m) {
    py::class_<SOMADenseNDArray>(m, "SOMADenseNDArray")

    .def_static("create", py::overload_cast<std::string_view, ArraySchema, std::map<std::string, std::string>>(&SOMADenseNDArray::create))
    .def_static("create", py::overload_cast<std::string_view, ArraySchema, std::shared_ptr<Context>>(&SOMADenseNDArray::create))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADenseNDArray::open))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::shared_ptr<Context>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADenseNDArray::open))

    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADenseNDArray::open))
    .def("close", &SOMADenseNDArray::close)
    .def("is_open", &SOMADenseNDArray::is_open)
    .def("type", &SOMADenseNDArray::type)
    .def("ctx", &SOMADenseNDArray::ctx)
    .def("is_sparse", &SOMADenseNDArray::is_sparse)
    .def("uri", &SOMADenseNDArray::uri)
    .def("schema", &SOMADenseNDArray::schema)
    .def("shape", &SOMADenseNDArray::shape)
    .def("ndim", &SOMADenseNDArray::ndim)
    .def("read_next", &SOMADenseNDArray::read_next)
    .def("write", &SOMADenseNDArray::write)
    .def("set_metadata", &SOMADenseNDArray::set_metadata)
    .def("delete_metadata", &SOMADenseNDArray::delete_metadata)
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMADenseNDArray::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMADenseNDArray::get_metadata))
    .def("has_metadata", &SOMADenseNDArray::has_metadata)
    .def("metadata_num", &SOMADenseNDArray::metadata_num);
}
}