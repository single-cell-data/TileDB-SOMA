/**
 * @file   soma_array.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <tiledbsoma/tiledbsoma>

#include "query_condition.cc"

#define DENUM(x) .value(#x, TILEDB_##x)

using namespace tiledbsoma;

namespace py = pybind11;
using namespace py::literals;

namespace tiledbsoma {

py::tuple get_enum(SOMAArray& sr, std::string attr_name){
    auto attr_to_enmrs = sr.get_attr_to_enum_mapping();
    if(attr_to_enmrs.count(attr_name) == 0)
        throw TileDBSOMAError("Given attribute does not have enumeration");

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
            throw TileDBSOMAError("Unsupported enumeration type.");
    }
}

bool get_enum_is_ordered(SOMAArray& sr, std::string attr_name){
    auto attr_to_enmrs = sr.get_attr_to_enum_mapping();
    if(attr_to_enmrs.count(attr_name) == 0)
        throw TileDBSOMAError("Given attribute does not have enumeration");
    return attr_to_enmrs.at(attr_name).ordered();
}

/**
 * @brief Convert ColumnBuffer to Arrow array.
 *
 * @param column_buffer ColumnBuffer
 * @return py::object Arrow array
 */
py::object to_array(std::shared_ptr<ColumnBuffer> column_buffer) {
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
py::object to_table(SOMAArray& sr, std::shared_ptr<ArrayBuffers> array_buffers) {
    auto pa = py::module::import("pyarrow");
    auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");
    auto pa_dict_from_arrays = pa.attr("DictionaryArray").attr("from_arrays");

    py::list names;
    py::list arrays;

    for (auto& name : array_buffers->names()) {
        auto column = array_buffers->at(name);
        names.append(name);

        if(sr.get_attr_to_enum_mapping().count(name) == 0){
            arrays.append(to_array(column));
        }else{
            arrays.append(pa_dict_from_arrays(
                to_array(column),
                get_enum(sr, name),
                py::none(),
                get_enum_is_ordered(sr, name)));
        }
    }

    auto pa_table = pa_table_from_arrays(arrays, names);

    return pa_table;
}

/**
 * @brief pybind11 bindings
 *
 */
void init_soma_array(py::module &m) {
    tiledbpy::init_query_condition(m);

    m.doc() = "SOMA acceleration library";

    m.def("version", []() { return tiledbsoma::version::as_string(); });

    m.def(
        "config_logging",
        [](const std::string& level, const std::string& logfile) {
            LOG_CONFIG(level, logfile);
        },
        "level"_a,
        "logfile"_a = "");

    m.def("info", &LOG_INFO, "message"_a = "");
    m.def("debug", &LOG_DEBUG, "message"_a = "");

    m.def(
        "tiledbsoma_stats_enable",
        []() { tiledbsoma::stats::enable(); },
        "Enable TileDB internal statistics. Lifecycle: experimental.");
    m.def(
        "tiledbsoma_stats_disable",
        []() { tiledbsoma::stats::disable(); },
        "Disable TileDB internal statistics. Lifecycle: experimental.");
    m.def(
        "tiledbsoma_stats_reset",
        []() { tiledbsoma::stats::reset(); },
        "Reset all TileDB internal statistics to 0. Lifecycle: experimental.");
    m.def(
        "tiledbsoma_stats_dump",
        []() {
            py::print(tiledbsoma::version::as_string());
            std::string stats = tiledbsoma::stats::dump();
            py::print(stats);
        },
        "Print TileDB internal statistics. Lifecycle: experimental.");

    py::class_<SOMAArray>(m, "SOMAArray")
        .def(
            py::init(
                [](std::string_view uri,
                   std::string_view name,
                   std::optional<std::vector<std::string>> column_names_in,
                   py::object py_query_condition,
                   std::string_view batch_size,
                   ResultOrder result_order,
                   std::map<std::string, std::string> platform_config,
                   std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                    // Handle optional args
                    std::vector<std::string> column_names;
                    if (column_names_in) {
                        column_names = *column_names_in;
                    }

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
                                init_pyqc(uri, column_names, platform_config, timestamp)
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
                                 .cast<tiledbpy::PyQueryCondition>()
                                 .ptr()
                                 .get();
                    }

                    // Release python GIL after we're done accessing python
                    // objects
                    py::gil_scoped_release release;

                    auto reader = SOMAArray::open(
                        OpenMode::read,
                        uri,
                        name,
                        platform_config,
                        column_names,
                        batch_size,
                        result_order,
                        timestamp);

                    // Set query condition if present
                    if (qc) {
                        reader->set_condition(*qc);
                    }

                    return reader;
                }),
            "uri"_a,
            py::kw_only(),
            "name"_a = "unnamed",
            "column_names"_a = py::none(),
            "query_condition"_a = py::none(),
            "batch_size"_a = "auto",
            "result_order"_a = ResultOrder::automatic,
            "platform_config"_a = py::dict(),
            "timestamp"_a = py::none())

        .def(
            "reset",
            [](SOMAArray& reader,
               std::optional<std::vector<std::string>> column_names_in,
               py::object py_query_condition,
               std::string_view batch_size,
               ResultOrder result_order) {
                // Handle optional args
                std::vector<std::string> column_names;
                if (column_names_in) {
                    column_names = *column_names_in;
                }

                // Handle query condition based on
                // TileDB-Py::PyQuery::set_attr_cond()
                QueryCondition* qc = nullptr;
                if (!py_query_condition.is(py::none())) {
                    py::object init_pyqc = py_query_condition.attr(
                        "init_query_condition");

                    try {
                        // Convert TileDB::Config to std::unordered map for pybind11 passing
                        std::unordered_map<std::string, std::string> cfg;
                        for (const auto& it : reader.ctx()->config()) {
                            cfg[it.first] = it.second;
                        }
                        // Column names will be updated with columns present in
                        // the query condition
                        auto new_column_names =
                            init_pyqc(reader.uri(), column_names, cfg, reader.timestamp())
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
                             .cast<tiledbpy::PyQueryCondition>()
                             .ptr()
                             .get();
                }

                // Release python GIL after we're done accessing python objects
                py::gil_scoped_release release;

                // Reset state of the existing SOMAArray object
                reader.reset(column_names, batch_size, result_order);

                // Set query condition if present
                if (qc) {
                    reader.set_condition(*qc);
                }
            },
            py::kw_only(),
            "column_names"_a = py::none(),
            "query_condition"_a = py::none(),
            "batch_size"_a = "auto",
            "result_order"_a = ResultOrder::automatic)

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
                        throw TileDBSOMAError(fmt::format(
                            "[pytiledbsoma] set_dim_points: type={} not "
                            "supported",
                            arrow_schema.format));
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
            "set_dim_points_float64",
            static_cast<void (SOMAArray::*)(
                const std::string&, const std::vector<double>&)>(
                &SOMAArray::set_dim_points))

        .def(
            "set_dim_points_float32",
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
            "set_dim_ranges_float64",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<double, double>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "set_dim_ranges_float32",
            static_cast<void (SOMAArray::*)(
                const std::string&,
                const std::vector<std::pair<float, float>>&)>(
                &SOMAArray::set_dim_ranges))

        .def(
            "submit",
            &SOMAArray::submit,
            py::call_guard<py::gil_scoped_release>())

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
                    return to_table(reader, *buffers);
                }

                // No data was read, the query is complete, return nullopt
                return std::nullopt;
            })

        .def("nnz", &SOMAArray::nnz, py::call_guard<py::gil_scoped_release>())

        .def_property_readonly("shape", &SOMAArray::shape)

        .def("get_enum", get_enum)

        .def("get_enum_is_ordered", get_enum_is_ordered)

        .def("get_enum_label_on_attr", &SOMAArray::get_enum_label_on_attr);
}
}  // namespace tiledbsoma