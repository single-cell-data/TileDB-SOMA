#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

using namespace tiledbsoma;

namespace py = pybind11;

namespace tiledbsoma {
void init_soma_experiment(py::module &m) {
    py::class_<SOMAExperiment>(m, "SOMAExperiment")

    .def_static("create", py::overload_cast<std::string_view, ArraySchema, std::map<std::string, std::string>>(&SOMAExperiment::create))
    .def_static("create", py::overload_cast<std::string_view, ArraySchema, std::shared_ptr<Context>>(&SOMAExperiment::create))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMAExperiment::open))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::shared_ptr<Context>, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMAExperiment::open))

    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMAExperiment::open))
    .def("close", &SOMAExperiment::close)
    .def("is_open", &SOMAExperiment::is_open)
    .def("type", &SOMAExperiment::type)
    .def("uri", &SOMAExperiment::uri)
    .def("ctx", &SOMAExperiment::ctx)
    .def("set", &SOMAExperiment::set)
    .def("get", &SOMAExperiment::get)
    .def("has", &SOMAExperiment::has)
    .def("count", &SOMAExperiment::count)
    .def("del", &SOMAExperiment::del)
    .def("member_to_uri_mapping", &SOMAExperiment::member_to_uri_mapping)
    .def("add_new_collection", &SOMAExperiment::add_new_collection)
    .def("add_new_experiment", &SOMAExperiment::add_new_experiment)
    .def("add_new_measurement", &SOMAExperiment::add_new_measurement)
    .def("add_new_dataframe", &SOMAExperiment::add_new_dataframe)
    .def("add_new_dense_ndarray", &SOMAExperiment::add_new_dense_ndarray)
    .def("add_new_sparse_ndarray", &SOMAExperiment::add_new_sparse_ndarray)
    .def("set_metadata", &SOMAExperiment::set_metadata)
    .def("delete_metadata", &SOMAExperiment::delete_metadata)
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMAExperiment::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMAExperiment::get_metadata))
    .def("has_metadata", &SOMAExperiment::has_metadata)
    .def("metadata_num", &SOMAExperiment::metadata_num);
}
}