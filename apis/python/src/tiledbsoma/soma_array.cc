/**
 * @file   soma_array.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
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
 * This file defines the SOMAArray bindings.
 */

#include "common.h"

#define DENUM(x) .value(#x, TILEDB_##x)
namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

py::tuple get_enum(SOMAArray& sr, std::string attr_name){
    auto attr_to_enmrs = sr.get_attr_to_enum_mapping();
    if(attr_to_enmrs.count(attr_name) == 0)
        TPY_ERROR_LOC("Given attribute does not have enumeration");

    Enumeration enmr(attr_to_enmrs.at(attr_name));

    switch (enmr.type()) {
        case TILEDB_UINT8:
            return py::tuple(py::cast(enmr.as_vector<uint8_t>()));
        case TILEDB_INT8:
            return py::tuple(py::cast(enmr.as_vector<int8_t>()));
        case TILEDB_UINT16:
            return py::tuple(py::cast(enmr.as_vector<uint16_t>()));
        case TILEDB_INT16:
            return py::tuple(py::cast(enmr.as_vector<int16_t>()));
        case TILEDB_UINT32:
            return py::tuple(py::cast(enmr.as_vector<uint32_t>()));
        case TILEDB_INT32:
            return py::tuple(py::cast(enmr.as_vector<int32_t>()));
        case TILEDB_UINT64:
            return py::tuple(py::cast(enmr.as_vector<uint64_t>()));
        case TILEDB_INT64:
            return py::tuple(py::cast(enmr.as_vector<int64_t>()));
        case TILEDB_FLOAT32:
            return py::tuple(py::cast(enmr.as_vector<float>()));
        case TILEDB_FLOAT64:
            return py::tuple(py::cast(enmr.as_vector<double>()));
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            return py::tuple(py::cast(enmr.as_vector<std::string>()));
        case TILEDB_BOOL:
            return py::tuple(py::cast(enmr.as_vector<bool>()));
        default:
            TPY_ERROR_LOC("Unsupported enumeration type.");
    }
}

bool get_enum_is_ordered(SOMAArray& sr, std::string attr_name){
    auto attr_to_enmrs = sr.get_attr_to_enum_mapping();
    if(attr_to_enmrs.count(attr_name) == 0)
        TPY_ERROR_LOC("Given attribute does not have enumeration");
    return attr_to_enmrs.at(attr_name).ordered();
}

void load_soma_array(py::module &m) {
    py::class_<SOMAArray>(m, "SOMAArray", "SOMAObject")
        .def(
            py::init(
                [](std::string_view uri,
                   std::string_view name,
                   std::optional<std::vector<std::string>> column_names_in,
                   std::string_view batch_size,
                   ResultOrder result_order,
                   std::map<std::string, std::string> platform_config,
                   std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                    // Handle optional args
                    std::vector<std::string> column_names;
                    if (column_names_in) {
                        column_names = *column_names_in;
                    }

                    return SOMAArray::open(
                        OpenMode::read,
                        uri,
                        name,
                        platform_config,
                        column_names,
                        batch_size,
                        result_order,
                        timestamp);
                }),
            "uri"_a,
            py::kw_only(),
            "name"_a = "unnamed",
            "column_names"_a = py::none(),
            "batch_size"_a = "auto",
            "result_order"_a = ResultOrder::automatic,
            "platform_config"_a = py::dict(),
            "timestamp"_a = py::none())

        .def(
            "set_condition", 
            [](SOMAArray& reader, 
               py::object py_query_condition,
               py::object py_schema){
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
                               init_pyqc(py_schema, column_names)
                                   .cast<std::vector<std::string>>();   
                           // Update the column_names list if it was not empty,
                           // otherwise continue selecting all columns with an
                           // empty column_names list
                           if (!column_names.empty()) {
                               column_names = new_column_names;
                           }
                       } catch (const std::exception& e) {
                           TPY_ERROR_LOC(e.what());
                       }   
                       qc = py_query_condition.attr("c_obj")
                                   .cast<PyQueryCondition>()
                                   .ptr()
                                   .get();
                   }   
                   reader.reset(column_names);

                    // Release python GIL after we're done accessing python
                   // objects
                   py::gil_scoped_release release;   
                   // Set query condition if present
                   if (qc) {
                       reader.set_condition(*qc);
                   }
                }, 
            "py_query_condition"_a,
            "py_schema"_a)

        .def(
            "reset",
            [](SOMAArray& reader,
               std::optional<std::vector<std::string>> column_names_in,
               std::string_view batch_size,
               ResultOrder result_order) {
                // Handle optional args
                std::vector<std::string> column_names;
                if (column_names_in) {
                    column_names = *column_names_in;
                }

                // Reset state of the existing SOMAArray object
                reader.reset(column_names, batch_size, result_order);
            },
            py::kw_only(),
            "column_names"_a = py::none(),
            "batch_size"_a = "auto",
            "result_order"_a = ResultOrder::automatic)
        
        .def("reopen", py::overload_cast<OpenMode, std::optional<std::pair<uint64_t, uint64_t>>>(&SOMAArray::open))
        .def("close", &SOMAArray::close)
        .def_property_readonly("closed", [](SOMAArray& reader) -> bool { 
            return not reader.is_open();
        })
        .def_property_readonly("mode", [](SOMAArray& reader){
            return reader.mode() == OpenMode::read ? "r" : "w";
        })
        .def_property_readonly("schema", [](SOMAArray& reader) -> py::object {
            auto pa = py::module::import("pyarrow");
            auto pa_schema_import = pa.attr("Schema").attr("_import_from_c");
            return pa_schema_import(py::capsule(reader.arrow_schema().get()));
        })
        .def("config", [](SOMAArray& reader) -> py::dict {
            return py::cast(reader.config());
        })

        // After this are short functions expected to be invoked when the coords
        // are Python list/tuple, or NumPy arrays.  Arrow arrays are in this
        // long if-else-if function.
        .def(
            "set_dim_points_arrow",
            [](SOMAArray& reader,
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
                        TPY_ERROR_LOC(
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

        // The following short functions are expected to be invoked when the
        // coords are Python list/tuple, or NumPy arrays.  Arrow arrays are in
        // the long if-else-if function above.
        //
        // Binding overloaded methods to templated member functions requires
        // more effort, see:
        // https://pybind11.readthedocs.io/en/stable/classes.html#overloaded-methods

        // In an initial version of this file we had `set_dim_ranges` relying
        // solely on type-overloading. This worked since we supported only int
        // and string indices. In a subsequent version we are now supporting
        // various NumPy/PyArrow types including float32, float64, int8, uint16,
        // etc. It is an unfortunate fact that pybind11 does _not_ successfully
        // disambiguate between float32 and float64, or between int8 and int64,
        // etc. given that we ask it to disambiguate using not just types but
        // std::vector of types or std::vector of std::pair of types.
        // Experiments have shown that when both float32 and float64 are
        // implemented with overloaded names to be differentiated solely by
        // type, pybind11 uses the _first found_. Therefore it is necessary for
        // us to no longer use common overloaded names.

        .def(
            "set_dim_points_string_or_bytes",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<std::string>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_double",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<double>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_float",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<float>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_int64",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<int64_t>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_int32",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<int32_t>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_int16",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<int16_t>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_int8",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<int8_t>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_uint64",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<uint64_t>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_uint32",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<uint32_t>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_uint16",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<uint16_t>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_uint8",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<uint8_t>&)>(
                &SOMAArray::set_dim_points))

        // In an initial version of this file we had `set_dim_ranges` relying
        // solely on type-overloading. This worked since we supported only int
        // and string indices. In a subsequent version we are now supporting
        // various NumPy/PyArrow types including float32, float64, int8, uint16,
        // etc. It is an unfortunate fact that pybind11 does _not_ successfully
        // disambiguate between float32 and float64, or between int8 and int64,
        // etc. given that we ask it to disambiguate using not just types but
        // std::vector of types or std::vector of std::pair of types.
        // Experiments have shown that when both float32 and float64 are
        // implemented with overloaded names to be differentiated solely by
        // type, pybind11 uses the _first found_. Therefore it is necessary for
        // us to no longer use common overloaded names.

        .def(
            "set_dim_ranges_string_or_bytes",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<std::string, std::string>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_int64",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<int64_t, int64_t>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_int32",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<int32_t, int32_t>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_int16",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<int16_t, int16_t>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_int8",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<int8_t, int8_t>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_uint64",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<uint64_t, uint64_t>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_uint32",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<uint32_t, uint32_t>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_uint16",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<uint16_t, uint16_t>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_uint8",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<uint8_t, uint8_t>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_double",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<double, double>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_float",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<float, float>>&)>(
                &SOMAArray::set_dim_ranges))

        .def("results_complete", &SOMAArray::results_complete)

        .def(
            "read_next",
            [](SOMAArray& reader) -> std::optional<py::object> {
                // Release python GIL before reading data
                py::gil_scoped_release release;

                // Try to read more data
                auto buffers = reader.read_next();

                // If more data was read, convert it to an arrow table and
                // return
                if (buffers.has_value()) {
                    // Acquire python GIL before accessing python objects
                    py::gil_scoped_acquire acquire;
                    return to_table(*buffers);
                }

                // No data was read, the query is complete, return nullopt
                return std::nullopt;
            })

        .def("nnz", &SOMAArray::nnz, py::call_guard<py::gil_scoped_release>())

        .def_property_readonly("shape", &SOMAArray::shape)

        .def_property_readonly("uri", &SOMAArray::uri)

        .def_property_readonly("column_names", &SOMAArray::column_names)

        .def_property_readonly("result_order", &SOMAArray::result_order)

        .def("get_enum", get_enum)

        .def("get_enum_is_ordered", get_enum_is_ordered)

        .def("get_enum_label_on_attr", &SOMAArray::get_enum_label_on_attr)

        .def_property_readonly("timestamp", [](SOMAArray& reader) -> py::object {
            if(!reader.timestamp().has_value())
                return py::none();
            return py::cast(reader.timestamp()->second);
        })

        .def("non_empty_domain", [](SOMAArray& reader, std::string name, py::dtype dtype){
            switch (np_to_tdb_dtype(dtype)) {
            case TILEDB_UINT64:
                return py::cast(reader.non_empty_domain<uint64_t>(name));
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
            case TILEDB_INT64:
                return py::cast(reader.non_empty_domain<int64_t>(name));
            case TILEDB_UINT32:
                return py::cast(reader.non_empty_domain<uint32_t>(name));
            case TILEDB_INT32:
                return py::cast(reader.non_empty_domain<int32_t>(name));
            case TILEDB_UINT16:
                return py::cast(reader.non_empty_domain<uint16_t>(name));
            case TILEDB_INT16:
                return py::cast(reader.non_empty_domain<int16_t>(name));
            case TILEDB_UINT8:
                return py::cast(reader.non_empty_domain<uint8_t>(name));
            case TILEDB_INT8:
                return py::cast(reader.non_empty_domain<int8_t>(name));
            case TILEDB_FLOAT64:
                return py::cast(reader.non_empty_domain<double>(name));
            case TILEDB_FLOAT32:
                return py::cast(reader.non_empty_domain<float>(name));
            case TILEDB_STRING_UTF8:
            case TILEDB_STRING_ASCII: 
                return py::cast(reader.non_empty_domain_var(name));
            default:
                throw TileDBSOMAError("Unsupported dtype for nonempty domain.");
            }
        })
        .def("domain", [](SOMAArray& reader, std::string name, py::dtype dtype) {
            switch (np_to_tdb_dtype(dtype)) {
            case TILEDB_UINT64:
                return py::cast(reader.domain<uint64_t>(name));
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
            case TILEDB_INT64:
                return py::cast(reader.domain<int64_t>(name));
            case TILEDB_UINT32:
                return py::cast(reader.domain<uint32_t>(name));
            case TILEDB_INT32:
                return py::cast(reader.domain<int32_t>(name));
            case TILEDB_UINT16:
                return py::cast(reader.domain<uint16_t>(name));
            case TILEDB_INT16:
                return py::cast(reader.domain<int16_t>(name));
            case TILEDB_UINT8:
                return py::cast(reader.domain<uint8_t>(name));
            case TILEDB_INT8:
                return py::cast(reader.domain<int8_t>(name));
            case TILEDB_FLOAT64:
                return py::cast(reader.domain<double>(name));
            case TILEDB_FLOAT32:
                return py::cast(reader.domain<float>(name));
            case TILEDB_STRING_UTF8:
            case TILEDB_STRING_ASCII: {
                std::pair<std::string, std::string> str_domain;
                return py::cast(std::make_pair("", ""));
            }
            default:
                throw TileDBSOMAError("Unsupported dtype for Dimension's domain");
            }
        })
        
        .def("set_metadata", &SOMAArray::set_metadata)
        .def("delete_metadata", &SOMAArray::delete_metadata)
        .def("get_metadata", 
            py::overload_cast<const std::string&>(&SOMAArray::get_metadata))
        .def_property_readonly("meta", [](SOMAArray&soma_dataframe) -> py::dict {
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
        })
        .def("has_metadata", &SOMAArray::has_metadata)
        .def("metadata_num", &SOMAArray::metadata_num);
}
}  // namespace tiledbsoma