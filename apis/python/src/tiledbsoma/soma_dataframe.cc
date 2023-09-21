#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

// #include <arrow/python/pyarrow.h>


using namespace tiledbsoma;

namespace py = pybind11;

// static std::unique_ptr<SOMADataFrame> create(
//         std::string_view uri, 
//         py::object pyschema, 
//         std::vector<std::string> index_column_names, 
//         std::map<std::string, std::string> platform_config, 
//         py::object pydomain){
//     auto schema_result = arrow::py::unwrap_schema(pyschema.ptr());
//     if(!schema_result.ok())
//         throw TileDBSOMAError(schema_result.status().ToString());
    
//     auto domain_result = arrow::py::unwrap_array(pydomain.ptr());
//     if(!domain_result.ok())
//         throw TileDBSOMAError(domain_result.status().ToString());

//     auto schema = schema_result.ValueOrDie();
//     auto domain = domain_result.ValueOrDie();

//     // struct ArrowSchema cschema {}; 
//     // SchemaExportGuard schema_guard(&cschema);
//     // auto status = ExportSchema(schema, &cschema);

//     // ExportArray(domain, cdomain);

//     return SOMADataFrame::create(
//         uri, cschema, index_column_names, platform_config, std::nullopt //TODO domain
//     );
// }

static std::unique_ptr<SOMADataFrame> create(
        std::string_view uri, 
        uintptr_t schema_ptr, 
        std::vector<std::string> index_column_names, 
        std::map<std::string, std::string> platform_config, 
        py::object pydomain){
    // auto schema_result = arrow::py::unwrap_schema(pyschema.ptr());
    // if(!schema_result.ok())
    //     throw TileDBSOMAError(schema_result.status().ToString());
    
    // auto domain_result = arrow::py::unwrap_array(pydomain.ptr());
    // if(!domain_result.ok())
    //     throw TileDBSOMAError(domain_result.status().ToString());

    // auto schema = schema_result.ValueOrDie();
    // auto domain = domain_result.ValueOrDie();

    // struct ArrowSchema cschema {}; 
    // SchemaExportGuard schema_guard(&cschema);
    // auto status = ExportSchema(schema, &cschema);

    // ExportArray(domain, cdomain);
    

    return SOMADataFrame::create(
        uri, (ArrowSchema*)schema_ptr, index_column_names, platform_config, std::nullopt //TODO domain
    );
}


namespace tiledbsoma {
void init_soma_dataframe(py::module &m) {
    // arrow::py::import_pyarrow();

    py::class_<SOMADataFrame>(m, "SOMADataFrame")

    .def_static("create", py::overload_cast<std::string_view,
        uintptr_t, // pointer to arrow schema
        std::vector<std::string>,
        std::map<std::string, std::string>,
        py::object>(create))

    // .def_static("create", py::overload_cast<std::string_view,
    //     ArrowSchema*,
    //     std::vector<std::string>,
    //     std::map<std::string, std::string>,
    //     std::optional<std::shared_ptr<ArrowArray>>>(&SOMADataFrame::create))
    // .def_static("create", py::overload_cast<std::string_view, py::object, std::vector<std::string>, std::map<std::string, std::string>, py::object>(create))
    // .def_static("create", py::overload_cast<std::string_view, py::capsule, std::shared_ptr<Context>>(&SOMADataFrame::create))
    // .def_static("create", py::overload_cast<std::string_view, struct ArrowSchema*, std::vector<std::string>, std::map<std::string, std::string>, std::optional<struct ArrowArray*>>(&SOMADataFrame::create))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::shared_ptr<Context>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))

    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))
    .def("close", &SOMADataFrame::close)
    .def("is_open", &SOMADataFrame::is_open)
    .def("type", &SOMADataFrame::type)
    .def("uri", &SOMADataFrame::uri)
    .def("ctx", &SOMADataFrame::ctx)
    .def("schema", &SOMADataFrame::schema)
    .def("index_column_names", &SOMADataFrame::index_column_names)
    .def("read_next", &SOMADataFrame::read_next)
    .def("write", &SOMADataFrame::write)
    .def("set_metadata", &SOMADataFrame::set_metadata)
    .def("delete_metadata", &SOMADataFrame::delete_metadata)
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMADataFrame::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMADataFrame::get_metadata))
    .def("has_metadata", &SOMADataFrame::has_metadata)
    .def("metadata_num", &SOMADataFrame::metadata_num);
}
}