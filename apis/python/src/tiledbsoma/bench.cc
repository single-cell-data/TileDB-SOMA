
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

#include "common.h"

namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_bench(py::module& m) {
    py::class_<MemoryBench>(m, "MemoryBench")
    .def_static(
        "alloc_vector",
        [](uint64_t size) {
            auto [pa_array, pa_schema] = MemoryBench::allocate_vector(size);

            std::cout << "Table constructed from C++" << std::endl;

            auto pa = py::module::import("pyarrow");
            auto pa_array_import = pa.attr("Array").attr("_import_from_c");
            auto pa_dtype_import = pa.attr("DataType").attr("_import_from_c");
            auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");

            py::list array_list;
            py::list field_list;

            std::string name(pa_schema->name);
            auto nullable = (pa_schema->flags & ARROW_FLAG_NULLABLE) != 0;
            auto dtype = pa_dtype_import(py::capsule(pa_schema.get()));
            array_list.append(pa_array_import(py::capsule(pa_array.get()), dtype));
            field_list.append(pa.attr("field")(std::string(name), dtype, nullable));
            
            return pa_table_from_arrays(array_list, "schema"_a = pa.attr("schema")(field_list));
        },
        py::kw_only(),
        "size"_a)
    .def_static(
        "alloc_pointer",
        [](uint64_t size) {
            auto [pa_array, pa_schema] = MemoryBench::allocate_pointer(size);

            std::cout << "Table constructed from C++" << std::endl;

            auto pa = py::module::import("pyarrow");
            auto pa_array_import = pa.attr("Array").attr("_import_from_c");
            auto pa_dtype_import = pa.attr("DataType").attr("_import_from_c");
            auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");

            py::list array_list;
            py::list field_list;

            std::string name(pa_schema->name);
            auto nullable = (pa_schema->flags & ARROW_FLAG_NULLABLE) != 0;
            auto dtype = pa_dtype_import(py::capsule(pa_schema.get()));
            array_list.append(pa_array_import(py::capsule(pa_array.get()), dtype));
            field_list.append(pa.attr("field")(std::string(name), dtype, nullable));
            
            return pa_table_from_arrays(array_list, "schema"_a = pa.attr("schema")(field_list));
        },
        py::kw_only(),
        "size"_a)
    .def_static(
        "alloc_column_buffer",
        [](uint64_t size) {
            auto [pa_array, pa_schema] = MemoryBench::allocate_column_buffer(size);

            std::cout << "Table constructed from C++" << std::endl;

            auto pa = py::module::import("pyarrow");
            auto pa_array_import = pa.attr("Array").attr("_import_from_c");
            auto pa_dtype_import = pa.attr("DataType").attr("_import_from_c");
            auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");

            py::list array_list;
            py::list field_list;

            std::string name(pa_schema->name);
            auto nullable = (pa_schema->flags & ARROW_FLAG_NULLABLE) != 0;
            auto dtype = pa_dtype_import(py::capsule(pa_schema.get()));
            array_list.append(pa_array_import(py::capsule(pa_array.get()), dtype));
            field_list.append(pa.attr("field")(std::string(name), dtype, nullable));
            
            return pa_table_from_arrays(array_list, "schema"_a = pa.attr("schema")(field_list));
        },
        py::kw_only(),
        "size"_a);
}
}