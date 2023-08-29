#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

using namespace tiledbsoma;

namespace py = pybind11;

namespace tiledbsoma {
void init_soma_collection(py::module &m) {
    py::class_<SOMACollection>(m, "SOMACollection")

    .def_static("create", py::overload_cast<std::string_view, std::map<std::string, std::string>>(&SOMACollection::create))
    .def_static("create", py::overload_cast<std::string_view, std::shared_ptr<Context>>(&SOMACollection::create))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMACollection::open))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::shared_ptr<Context>, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMACollection::open))

    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMACollection::open))
    .def("close", &SOMACollection::close)
    .def("is_open", &SOMACollection::is_open)
    .def("type", &SOMACollection::type)
    .def("uri", &SOMACollection::uri)
    .def("ctx", &SOMACollection::ctx)
    .def("set", &SOMACollection::set)
    .def("get", &SOMACollection::get)
    .def("has", &SOMACollection::has)
    .def("count", &SOMACollection::count)
    .def("del", &SOMACollection::del)
    .def("member_to_uri_mapping", &SOMACollection::member_to_uri_mapping)
    .def("add_new_collection", &SOMACollection::add_new_collection)
    .def("add_new_experiment", &SOMACollection::add_new_experiment)
    .def("add_new_measurement", &SOMACollection::add_new_measurement)
    .def("add_new_dataframe", &SOMACollection::add_new_dataframe)
    .def("add_new_dense_ndarray", &SOMACollection::add_new_dense_ndarray)
    .def("add_new_sparse_ndarray", &SOMACollection::add_new_sparse_ndarray)
    .def("set_metadata", &SOMACollection::set_metadata)
    .def("delete_metadata", &SOMACollection::delete_metadata)
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMACollection::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMACollection::get_metadata))
    .def("has_metadata", &SOMACollection::has_metadata)
    .def("metadata_num", &SOMACollection::metadata_num);
}
}