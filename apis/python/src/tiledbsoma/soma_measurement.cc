#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

using namespace tiledbsoma;

namespace py = pybind11;

namespace tiledbsoma {
void init_soma_measurement(py::module &m) {
    py::class_<SOMAMeasurement>(m, "SOMAMeasurement")

    .def_static("create", py::overload_cast<std::string_view, ArraySchema, std::map<std::string, std::string>>(&SOMAMeasurement::create))
    .def_static("create", py::overload_cast<std::string_view, ArraySchema, std::shared_ptr<Context>>(&SOMAMeasurement::create))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMAMeasurement::open))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::shared_ptr<Context>, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMAMeasurement::open))

    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMAMeasurement::open))
    .def("close", &SOMAMeasurement::close)
    .def("is_open", &SOMAMeasurement::is_open)
    .def("type", &SOMAMeasurement::type)
    .def("uri", &SOMAMeasurement::uri)
    .def("ctx", &SOMAMeasurement::ctx)
    .def("set", &SOMAMeasurement::set)
    .def("get", &SOMAMeasurement::get)
    .def("has", &SOMAMeasurement::has)
    .def("count", &SOMAMeasurement::count)
    .def("del", &SOMAMeasurement::del)
    .def("member_to_uri_mapping", &SOMAMeasurement::member_to_uri_mapping)
    .def("add_new_collection", &SOMAMeasurement::add_new_collection)
    .def("add_new_experiment", &SOMAMeasurement::add_new_experiment)
    .def("add_new_measurement", &SOMAMeasurement::add_new_measurement)
    .def("add_new_dataframe", &SOMAMeasurement::add_new_dataframe)
    .def("add_new_dense_ndarray", &SOMAMeasurement::add_new_dense_ndarray)
    .def("add_new_sparse_ndarray", &SOMAMeasurement::add_new_sparse_ndarray)
    .def("set_metadata", &SOMAMeasurement::set_metadata)
    .def("delete_metadata", &SOMAMeasurement::delete_metadata)
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMAMeasurement::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMAMeasurement::get_metadata))
    .def("has_metadata", &SOMAMeasurement::has_metadata)
    .def("metadata_num", &SOMAMeasurement::metadata_num);
}
}