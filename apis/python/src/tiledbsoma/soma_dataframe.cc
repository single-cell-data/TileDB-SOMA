#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

using namespace tiledbsoma;
namespace py = pybind11;

namespace tiledbsoma {

static std::unique_ptr<SOMADataFrame> create(
    std::string_view uri, 
    uintptr_t schema_ptr, 
    std::vector<std::string> index_column_names, 
    std::map<std::string, std::string> platform_config, 
    uintptr_t domains_ptr,
    uintptr_t extents_ptr){
    return SOMADataFrame::create(
        uri, *((ArrowSchema*)schema_ptr), platform_config, index_column_names, 
        *((ArrowArray*)domains_ptr), *((ArrowArray*)extents_ptr));
}

static void write(SOMADataFrame& dataframe, uintptr_t schema_ptr, uintptr_t array_ptr){
    ArrowSchema arrow_schema = *((ArrowSchema*)schema_ptr);
    ArrowArray arrow_array = *((ArrowArray*)array_ptr);
    auto schema = dataframe.schema();
    auto array_buffer = std::make_shared<ArrayBuffers>();

    for(int64_t i = 0; i < arrow_array.n_children; ++i){
        auto child = arrow_schema.children[i];
        auto name = child->name;
        auto typeinfo = ArrowAdapter::arrow_type_to_tiledb(child);
        array_buffer->emplace(name, 
            ColumnBuffer::create(
                *schema, 
                name, 
                arrow_array.children[i]->buffers[1], 
                arrow_array.children[i]->length * typeinfo.elem_size)); 
    }
    dataframe.write(array_buffer);
}

void init_soma_dataframe(py::module &m) {
    py::class_<SOMADataFrame>(m, "SOMADataFrame")

    .def_static("create", py::overload_cast<std::string_view,
        uintptr_t,
        std::vector<std::string>,
        std::map<std::string, std::string>,
        uintptr_t,
        uintptr_t>(create))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::shared_ptr<Context>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))
    .def_static("exists", &SOMADataFrame::exists)

    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))
    .def("close", &SOMADataFrame::close)
    .def("is_open", &SOMADataFrame::is_open)
    .def("type", &SOMADataFrame::type)
    .def("uri", &SOMADataFrame::uri)
    .def("ctx", &SOMADataFrame::ctx)
    .def("schema", &SOMADataFrame::schema)
    .def("index_column_names", &SOMADataFrame::index_column_names)
    .def("count", &SOMADataFrame::count)
    .def("read_next", &SOMADataFrame::read_next)
    .def("write", write)
    .def("set_metadata", &SOMADataFrame::set_metadata)
    .def("delete_metadata", &SOMADataFrame::delete_metadata)
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMADataFrame::get_metadata))
    .def("get_metadata", py::overload_cast<>(&SOMADataFrame::get_metadata))
    .def("has_metadata", &SOMADataFrame::has_metadata)
    .def("metadata_num", &SOMADataFrame::metadata_num);
}
}