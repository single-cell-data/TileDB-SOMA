#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

using namespace tiledbsoma;

namespace py = pybind11;

namespace tiledbsoma {
void init_soma_sparse_ndarray(py::module &m) {
    py::class_<SOMASparseNDArray>(m, "SOMASparseNDArray")

    .def_static("create", py::overload_cast<std::string_view, ArraySchema, std::map<std::string, std::string>>(&SOMASparseNDArray::create))
    .def_static("create", py::overload_cast<std::string_view, ArraySchema, std::shared_ptr<Context>>(&SOMASparseNDArray::create))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMASparseNDArray::open))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::shared_ptr<Context>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMASparseNDArray::open))

    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMASparseNDArray::open))
    .def("close", &SOMASparseNDArray::close)
    .def("is_open", &SOMASparseNDArray::is_open)
    .def("type", &SOMASparseNDArray::type)
    .def("ctx", &SOMASparseNDArray::ctx)
    .def("is_sparse", &SOMASparseNDArray::is_sparse)
    .def("uri", &SOMASparseNDArray::uri)
    .def("schema", &SOMASparseNDArray::schema)
    .def("shape", &SOMASparseNDArray::shape)
    .def("ndim", &SOMASparseNDArray::ndim)
    .def("nnz", &SOMASparseNDArray::nnz)
    .def("read_next", &SOMASparseNDArray::read_next)
    .def("write", &SOMASparseNDArray::write)
    .def("set_metadata", &SOMASparseNDArray::set_metadata)
    .def("delete_metadata", &SOMASparseNDArray::delete_metadata)
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMASparseNDArray::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMASparseNDArray::get_metadata))
    .def("has_metadata", &SOMASparseNDArray::has_metadata)
    .def("metadata_num", &SOMASparseNDArray::metadata_num);
}
}