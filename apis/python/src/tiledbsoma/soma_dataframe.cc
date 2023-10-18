/**
 * @file   soma_dataframe.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file defines the SOMADataFrame bindings.
 */

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

#include "common.h"

namespace py = pybind11;
using namespace py::literals;

namespace tiledbsoma {

/**
 * @brief Convert ArrayBuffers to Arrow table.
 *
 * @param cbs ArrayBuffers
 * @return py::object
 */
py::object df_to_table(std::shared_ptr<ArrayBuffers> array_buffers) {
    auto pa = py::module::import("pyarrow");
    auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");
    auto pa_array_import = pa.attr("Array").attr("_import_from_c");

    py::list array_list;
    py::list names;

    for (auto& name : array_buffers->names()) {
        auto column = array_buffers->at(name);
        auto [pa_array, pa_schema] = ArrowAdapter::to_arrow(column);
        auto array = pa_array_import(py::capsule(pa_array.get()), 
                                     py::capsule(pa_schema.get()));
        array_list.append(array);
        names.append(name);
    }

    return pa_table_from_arrays(array_list, names);
}

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

std::unordered_map<tiledb_datatype_t, std::string> _tdb_to_np_name_dtype = {
    {TILEDB_INT32, "int32"},
    {TILEDB_INT64, "int64"},
    {TILEDB_FLOAT32, "float32"},
    {TILEDB_FLOAT64, "float64"},
    {TILEDB_INT8, "int8"},
    {TILEDB_UINT8, "uint8"},
    {TILEDB_INT16, "int16"},
    {TILEDB_UINT16, "uint16"},
    {TILEDB_UINT32, "uint32"},
    {TILEDB_UINT64, "uint64"},
    {TILEDB_STRING_ASCII, "S"},
    {TILEDB_STRING_UTF8, "U1"},
    {TILEDB_CHAR, "S1"},
    {TILEDB_DATETIME_YEAR, "M8[Y]"},
    {TILEDB_DATETIME_MONTH, "M8[M]"},
    {TILEDB_DATETIME_WEEK, "M8[W]"},
    {TILEDB_DATETIME_DAY, "M8[D]"},
    {TILEDB_DATETIME_HR, "M8[h]"},
    {TILEDB_DATETIME_MIN, "M8[m]"},
    {TILEDB_DATETIME_SEC, "M8[s]"},
    {TILEDB_DATETIME_MS, "M8[ms]"},
    {TILEDB_DATETIME_US, "M8[us]"},
    {TILEDB_DATETIME_NS, "M8[ns]"},
    {TILEDB_DATETIME_PS, "M8[ps]"},
    {TILEDB_DATETIME_FS, "M8[fs]"},
    {TILEDB_DATETIME_AS, "M8[as]"},
    {TILEDB_TIME_HR, "m8[h]"},
    {TILEDB_TIME_MIN, "m8[m]"},
    {TILEDB_TIME_SEC, "m8[s]"},
    {TILEDB_TIME_MS, "m8[ms]"},
    {TILEDB_TIME_US, "m8[us]"},
    {TILEDB_TIME_NS, "m8[ns]"},
    {TILEDB_TIME_PS, "m8[ps]"},
    {TILEDB_TIME_FS, "m8[fs]"},
    {TILEDB_TIME_AS, "m8[as]"},
    {TILEDB_BLOB, "byte"},
    {TILEDB_BOOL, "bool"},
};

py::dtype tdb_to_np_dtype(tiledb_datatype_t type, uint32_t cell_val_num) {
  if (type == TILEDB_CHAR || type == TILEDB_STRING_UTF8 ||
      type == TILEDB_STRING_ASCII) {
    std::string base_str = (type == TILEDB_STRING_UTF8) ? "|U" : "|S";
    if (cell_val_num < TILEDB_VAR_NUM)
      base_str += std::to_string(cell_val_num);
    return py::dtype(base_str);
  }

  if (cell_val_num == 1) {
    if (type == TILEDB_STRING_UTF16 || type == TILEDB_STRING_UTF32)
      throw std::invalid_argument("Unimplemented UTF16 or UTF32 string conversion!");
    if (type == TILEDB_STRING_UCS2 || type == TILEDB_STRING_UCS4)
      throw std::invalid_argument("Unimplemented UCS2 or UCS4 string conversion!");

    if (_tdb_to_np_name_dtype.count(type) == 1)
      return py::dtype(_tdb_to_np_name_dtype[type]);
  }

  if (cell_val_num == 2) {
    if (type == TILEDB_FLOAT32)
      return py::dtype("complex64");
    if (type == TILEDB_FLOAT64)
      return py::dtype("complex128");
  }

  if (cell_val_num == TILEDB_VAR_NUM)
    return tdb_to_np_dtype(type, 1);

  if (cell_val_num > 1) {
    py::dtype base_dtype = tdb_to_np_dtype(type, 1);
    py::tuple rec_elem = py::make_tuple("", base_dtype);
    py::list rec_list;
    for (size_t i = 0; i < cell_val_num; i++)
      rec_list.append(rec_elem);
    // note: we call the 'dtype' constructor b/c py::dtype does not accept
    // list
    auto np = py::module::import("numpy");
    auto np_dtype = np.attr("dtype");
    return np_dtype(rec_list);
  }

  throw std::invalid_argument("tiledb datatype not understood ('" +
        tiledb::impl::type_to_str(type) +
        "', cell_val_num: " + std::to_string(cell_val_num) + ")");
}

py::dict meta(SOMADataFrame &soma_dataframe) {
    py::dict results;
    
    for (auto const& [key, val] : soma_dataframe.get_metadata()){
        tiledb_datatype_t tdb_type = std::get<MetadataInfo::dtype>(val);
        uint32_t value_num = std::get<MetadataInfo::num>(val);
        const void *value = std::get<MetadataInfo::value>(val);

        if(tdb_type == TILEDB_STRING_UTF8){
            results[py::str(key)] = py::str(std::string((const char*)value, value_num));
        }else if(tdb_type == TILEDB_STRING_ASCII){
            results[py::str(key)] = py::bytes(std::string((const char*)value, value_num));
        }else{
            py::dtype value_type = tdb_to_np_dtype(tdb_type, 1);
            results[py::str(key)] = py::array(value_type, value_num, value);
        }
    }

    return results;
}

void init_soma_dataframe(py::module &m) {
    py::class_<SOMADataFrame>(m, "SOMADataFrame")

    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::map<std::string, std::string>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))
    .def_static("open", py::overload_cast<std::string_view, OpenMode, std::shared_ptr<Context>, std::vector<std::string>, ResultOrder, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))
    .def_static("exists", &SOMADataFrame::exists)

    .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMADataFrame::open))
    .def("close", &SOMADataFrame::close)
    .def_property_readonly("closed", [](SOMADataFrame& soma_df) -> bool { 
        return soma_df.is_open();
    })
    .def("reset", &SOMADataFrame::reset)

    .def("set_condition", 
        [](SOMADataFrame& reader, 
            py::object py_query_condition,
            py::object pa_schema){  
                auto column_names = reader.column_names();
                // Handle query condition based on
                // TileDB-Py::PyQuery::set_attr_cond()
                QueryCondition* qc = nullptr;
                if (!py_query_condition.is(py::none())) {
                py::object init_pyqc = py_query_condition.attr(
                    "init_query_condition");   
                try {
                    // Column names will be updated with columns present
                    // in the query condition
                    auto new_column_names =
                        init_pyqc(pa_schema, column_names)
                            .cast<std::vector<std::string>>();   
                    // Update the column_names list if it was not empty,
                    // otherwise continue selecting all columns with an
                    // empty column_names list
                    if (!column_names.empty()) {
                        column_names = new_column_names;
                    }
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }   
                qc = py_query_condition.attr("c_obj")
                    .cast<PyQueryCondition>()
                    .ptr()
                    .get();
                reader.reset(column_names);

                // Release python GIL after we're done accessing python
                // objects
                py::gil_scoped_release release;   
                // Set query condition if present
                if (qc) {
                    reader.set_condition(*qc);
                }
            }
        }, 
        "py_query_condition"_a,
        "py_schema"_a)

    .def("type", &SOMADataFrame::type)
    .def("uri", &SOMADataFrame::uri)
    .def_property_readonly("schema", [](SOMADataFrame& soma_df) -> py::object {
        auto pa = py::module::import("pyarrow");
        auto pa_schema_import = pa.attr("Schema").attr("_import_from_c");
        auto schema = soma_df.schema();
        auto attr_to_enum = soma_df.get_attr_to_enum_mapping();
        auto tdb_schema = 
            ArrowAdapter::tiledb_schema_to_arrow_schema(schema, attr_to_enum);
        return pa_schema_import(py::capsule(tdb_schema.get()));
    })
    .def_property_readonly("timestamp", &SOMADataFrame::timestamp)
    .def_property_readonly("index_column_names", &SOMADataFrame::index_column_names)
    .def("nonempty_domain", [](SOMADataFrame& soma_df, std::string name){
        switch (soma_df.schema()->domain().dimension(name).type()) {
        case TILEDB_UINT64: {
            return py::cast(soma_df.non_empty_domain<uint64_t>(name));
        }
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_INT64: {
            return py::cast(soma_df.non_empty_domain<int64_t>(name));
        }
        case TILEDB_UINT32: {
            return py::cast(soma_df.non_empty_domain<uint32_t>(name));
        }
        case TILEDB_INT32: {
            return py::cast(soma_df.non_empty_domain<int32_t>(name));
        }
        case TILEDB_UINT16: {
            return py::cast(soma_df.non_empty_domain<uint16_t>(name));
        }
        case TILEDB_INT16: {
            return py::cast(soma_df.non_empty_domain<int16_t>(name));
        }
        case TILEDB_UINT8: {
            return py::cast(soma_df.non_empty_domain<uint8_t>(name));
        }
        case TILEDB_INT8: {
            return py::cast(soma_df.non_empty_domain<int8_t>(name));
        }
        case TILEDB_FLOAT64: {
            return py::cast(soma_df.non_empty_domain<double>(name));
        }
        case TILEDB_FLOAT32: {
            return py::cast(soma_df.non_empty_domain<float>(name));
        }
        case TILEDB_STRING_UTF8:
        case TILEDB_STRING_ASCII: {
            return py::cast(soma_df.non_empty_domain_var(name));
        }
        default:
            throw TileDBSOMAError("Unsupported dtype for nonempty domain.");
        }
    })
    .def("domain", [](SOMADataFrame& soma_df, std::string name) -> py::tuple {
        auto dim = soma_df.schema()->domain().dimension(name);
        switch (dim.type()) {
        case TILEDB_UINT64: {
        auto dom = dim.domain<uint64_t>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_INT64: {
        auto dom = dim.domain<int64_t>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_UINT32: {
        auto dom = dim.domain<uint32_t>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_INT32: {
        auto dom = dim.domain<int32_t>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_UINT16: {
        auto dom = dim.domain<uint16_t>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_INT16: {
        auto dom = dim.domain<int16_t>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_UINT8: {
        auto dom = dim.domain<uint8_t>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_INT8: {
        auto dom = dim.domain<int8_t>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_FLOAT64: {
        auto dom = dim.domain<double>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_FLOAT32: {
        auto dom = dim.domain<float>();
        return py::make_tuple(dom.first, dom.second);
        }
        case TILEDB_STRING_ASCII: {
        return py::make_tuple("", "");
        }
        default:
        throw TileDBSOMAError("Unsupported dtype for Dimension's domain");
        }
    })
    .def_property_readonly("ndim", &SOMADataFrame::ndim)
    .def_property_readonly("count", &SOMADataFrame::count)
    .def("read_next", [](SOMADataFrame& dataframe){
        // Release GIL when reading data
        py::gil_scoped_release release;
        auto buffers = dataframe.read_next();
        py::gil_scoped_acquire acquire;

        return to_table(buffers);
    })
    
    .def("set_metadata", &SOMADataFrame::set_metadata)
    .def("delete_metadata", &SOMADataFrame::delete_metadata)
    .def("get_metadata", 
        py::overload_cast<const std::string&>(&SOMADataFrame::get_metadata))
    .def_property_readonly("meta", [](SOMADataFrame&soma_dataframe) -> py::dict {
        py::dict results;
            
            for (auto const& [key, val] : soma_dataframe.get_metadata()){
                tiledb_datatype_t tdb_type = std::get<MetadataInfo::dtype>(val);
                uint32_t value_num = std::get<MetadataInfo::num>(val);
                const void *value = std::get<MetadataInfo::value>(val);

                if(tdb_type == TILEDB_STRING_UTF8){
                    results[py::str(key)] = py::str(std::string((const char*)value, value_num));
                }else if(tdb_type == TILEDB_STRING_ASCII){
                    results[py::str(key)] = py::bytes(std::string((const char*)value, value_num));
                }else{
                    py::dtype value_type = tdb_to_np_dtype(tdb_type, 1);
                    results[py::str(key)] = py::array(value_type, value_num, value);
                }
            }

            return results;
        }
    )
    .def("has_metadata", &SOMADataFrame::has_metadata)
    .def("metadata_num", &SOMADataFrame::metadata_num)
    .def(
        "set_dim_points_arrow",
        [](SOMADataFrame& reader,
            const std::string& dim,
            py::object py_arrow_array,
            int partition_index,
            int partition_count) {
            // Create a list of array chunks
            py::list array_chunks;
            if (py::hasattr(py_arrow_array, "chunks")) {
                array_chunks = py_arrow_array.attr("chunks")
                                    .cast<py::list>();
            } else {
                array_chunks.append(py_arrow_array);
            }

            for (const pybind11::handle array : array_chunks) {
                ArrowSchema arrow_schema;
                ArrowArray arrow_array;
                uintptr_t arrow_schema_ptr = (uintptr_t)(&arrow_schema);
                uintptr_t arrow_array_ptr = (uintptr_t)(&arrow_array);

                // Call array._export_to_c to get arrow array and schema
                //
                // If ever a NumPy array gets in here, there will be an
                // exception like "AttributeError: 'numpy.ndarray' object
                // has no attribute '_export_to_c'".
                array.attr("_export_to_c")(
                    arrow_array_ptr, arrow_schema_ptr);

                auto coords = array.attr("tolist")();

                if (!strcmp(arrow_schema.format, "l")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<int64_t>>());
                } else if (!strcmp(arrow_schema.format, "i")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<int32_t>>());
                } else if (!strcmp(arrow_schema.format, "s")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<int16_t>>());
                } else if (!strcmp(arrow_schema.format, "c")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<int8_t>>());
                } else if (!strcmp(arrow_schema.format, "L")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<uint64_t>>());
                } else if (!strcmp(arrow_schema.format, "I")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<uint32_t>>());
                } else if (!strcmp(arrow_schema.format, "S")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<uint16_t>>());
                } else if (!strcmp(arrow_schema.format, "C")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<uint8_t>>());
                } else if (!strcmp(arrow_schema.format, "f")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<float>>());
                } else if (!strcmp(arrow_schema.format, "g")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<double>>());
                } else if (
                    !strcmp(arrow_schema.format, "u") ||
                    !strcmp(arrow_schema.format, "z")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<std::string>>());
                } else if (
                    !strcmp(arrow_schema.format, "tss:") ||
                    !strcmp(arrow_schema.format, "tsm:") ||
                    !strcmp(arrow_schema.format, "tsu:") ||
                    !strcmp(arrow_schema.format, "tsn:")) {
                    // convert the Arrow Array to int64
                    auto pa = py::module::import("pyarrow");
                    coords = array.attr("cast")(pa.attr("int64")()).attr("tolist")();
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<int64_t>>());
                } else if (
                    !strcmp(arrow_schema.format, "U") ||
                    !strcmp(arrow_schema.format, "Z")) {
                    reader.set_dim_points(
                        dim, coords.cast<std::vector<std::string>>());
                } else {
                    throw TileDBSOMAError(
                        "[pytiledbsoma] set_dim_points: type={} not "
                        "supported" + 
                        std::string(arrow_schema.format));
                }

                // Release arrow schema
                arrow_schema.release(&arrow_schema);
            }
        },
        "dim"_a,
        "py_arrow_array"_a,
        "partition_index"_a = 0,
        "partition_count"_a = 1)

    .def(
        "set_dim_points_string_or_bytes",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<std::string>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_float64",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<double>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_float32",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<float>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_int64",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<int64_t>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_int32",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<int32_t>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_int16",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<int16_t>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_int8",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<int8_t>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_uint64",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<uint64_t>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_uint32",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<uint32_t>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_uint16",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<uint16_t>&)>(
            &SOMADataFrame::set_dim_points))
    .def(
        "set_dim_points_uint8",
        static_cast<void (SOMADataFrame::*)(
            const std::string&, const std::vector<uint8_t>&)>(
            &SOMADataFrame::set_dim_points))
            
    .def(
        "set_dim_ranges_string_or_bytes",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<std::string, std::string>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_int64",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<int64_t, int64_t>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_int32",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<int32_t, int32_t>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_int16",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<int16_t, int16_t>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_int8",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<int8_t, int8_t>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_uint64",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<uint64_t, uint64_t>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_uint32",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<uint32_t, uint32_t>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_uint16",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<uint16_t, uint16_t>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_uint8",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<uint8_t, uint8_t>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_float64",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<double, double>>&)>(
            &SOMADataFrame::set_dim_ranges))

    .def(
        "set_dim_ranges_float32",
        static_cast<void (SOMADataFrame::*)(
            const std::string&,
            const std::vector<std::pair<float, float>>&)>(
            &SOMADataFrame::set_dim_ranges));
}
}