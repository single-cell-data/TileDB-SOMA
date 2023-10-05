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
    uintptr_t ptr_dom_and_exts)
{
    return SOMADataFrame::create(
        uri, *((ArrowSchema*)schema_ptr), platform_config, index_column_names, *((ArrowArray*)ptr_dom_and_exts));
}

/**
 * @brief Convert ColumnBuffer to Arrow array.
 *
 * @param column_buffer ColumnBuffer
 * @return py::object Arrow array
 */
py::object df_to_array(std::shared_ptr<ColumnBuffer> column_buffer) {
    auto pa = py::module::import("pyarrow");
    auto pa_array_import = pa.attr("Array").attr("_import_from_c");

    auto [array, schema] = ArrowAdapter::to_arrow(column_buffer);
    return pa_array_import(py::capsule(array.get()), py::capsule(schema.get()));
}

/**
 * @brief Convert ArrayBuffers to Arrow table.
 *
 * @param cbs ArrayBuffers
 * @return py::object
 */
py::object df_to_table(std::shared_ptr<ArrayBuffers> array_buffers) {
    auto pa = py::module::import("pyarrow");
    auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");
    auto pa_dict_from_arrays = pa.attr("DictionaryArray").attr("from_arrays");

    py::list names;
    py::list arrays;

    for (auto& name : array_buffers->names()) {
        auto column = array_buffers->at(name);
        names.append(name);
        arrays.append(df_to_array(column));
    }

    auto pa_table = pa_table_from_arrays(arrays, names);

    return pa_table;
}

static void write(SOMADataFrame& dataframe, 
    uintptr_t schema_ptr, 
    uintptr_t array_ptr)
{
    ArrowSchema arrow_schema = *((ArrowSchema*)schema_ptr);
    ArrowArray arrow_array = *((ArrowArray*)array_ptr);
    auto array_buffer = std::make_shared<ArrayBuffers>();

    for(int64_t i = 0; i < arrow_array.n_children; ++i){
        auto schema_child = arrow_schema.children[i];
        auto array_child = arrow_array.children[i]; // TODO DEAL WITH STRINGS
        auto schema = dataframe.schema();
        auto typeinfo = ArrowAdapter::arrow_type_to_tiledb(schema_child);

        std::vector<int64_t> data;
        data.assign(
            (const int64_t*)array_child->buffers[1],
            (const int64_t*)array_child->buffers[1] + 5
        );
        std::cout << "DF " << schema_child->name << " " << schema_child->format << " " << array_child->length << std::endl;
        for(auto d : data)
            std::cout << d << " ";
        std::cout << std::endl;

        array_buffer->emplace(schema_child->name, 
            ColumnBuffer::create(
                *schema, 
                schema_child->name, 
                array_child->buffers[1], 
                array_child->length,
                typeinfo.elem_size)); 
    }
    dataframe.write(array_buffer);
}

// TODO RETURN POINTERS, CONSTURCT PYARROW TABLE IN PYTHON CODE
static std::optional<py::object> read_next(SOMADataFrame& dataframe){
    // Release python GIL before reading data
    py::gil_scoped_release release;

    // Try to read more data
    auto buffers = dataframe.read_next();

    // If more data was read, convert it to an arrow table and return
    if (buffers.has_value()) {
        // Acquire python GIL before accessing python objects
        py::gil_scoped_acquire acquire;
        return df_to_table(*buffers);
    }

    // No data was read, the query is complete, return nullopt
    return std::nullopt;
}

void init_soma_dataframe(py::module &m) {
    py::class_<SOMADataFrame>(m, "SOMADataFrame")

    .def_static("create", py::overload_cast<std::string_view,
        uintptr_t,
        std::vector<std::string>,
        std::map<std::string, std::string>,
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
    .def("read_next", read_next)
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