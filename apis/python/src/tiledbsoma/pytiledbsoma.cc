/**
 * @file   pytiledbsoma.cc
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
 * This file defines the a pybind11 api into SOMA C++ library.
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
py::object to_table(std::shared_ptr<ArrayBuffers> array_buffers) {
    auto pa = py::module::import("pyarrow");
    auto pa_table_from_arrays = pa.attr("Table").attr("from_arrays");

    py::list names;
    py::list arrays;

    for (auto& name : array_buffers->names()) {
        auto column = array_buffers->at(name);
        names.append(name);
        arrays.append(to_array(column));
    }

    return pa_table_from_arrays(arrays, names);
}

/**
 * @brief pybind11 bindings
 *
 */
PYBIND11_MODULE(pytiledbsoma, m) {
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
                   py::object py_schema,
                   std::string_view batch_size,
                   std::string_view result_order,
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
                                init_pyqc(py_schema, column_names)
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
                        TILEDB_READ,
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
            "schema"_a = py::none(),
            "batch_size"_a = "auto",
            "result_order"_a = "auto",
            "platform_config"_a = py::dict(),
            "timestamp"_a = py::none())

        .def(
            "reset",
            [](SOMAArray& reader,
               std::optional<std::vector<std::string>> column_names_in,
               py::object py_query_condition,
               py::object py_schema,
               std::string_view batch_size,
               std::string_view result_order) {
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
                        // Column names will be updated with columns present in
                        // the query condition
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
            "schema"_a = py::none(),
            "batch_size"_a = "auto",
            "result_order"_a = "auto")

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

                    int data_index = arrow_array.n_buffers - 1;

                    if (!strcmp(arrow_schema.format, "l")) {
                        tcb::span<int64_t> data{
                            (int64_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "i")) {
                        tcb::span<int32_t> data{
                            (int32_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "s")) {
                        tcb::span<int16_t> data{
                            (int16_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "c")) {
                        tcb::span<int8_t> data{
                            (int8_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "L")) {
                        tcb::span<uint64_t> data{
                            (uint64_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "I")) {
                        tcb::span<uint32_t> data{
                            (uint32_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "S")) {
                        tcb::span<uint16_t> data{
                            (uint16_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "C")) {
                        tcb::span<uint8_t> data{
                            (uint8_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "f")) {
                        tcb::span<float> data{
                            (float*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (!strcmp(arrow_schema.format, "g")) {
                        tcb::span<double> data{
                            (double*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                    } else if (
                        !strcmp(arrow_schema.format, "u") ||
                        !strcmp(arrow_schema.format, "z")) {
                        // TODO: partitioning is not supported for string/bytes
                        // dims
                        const char* data = (const char*)(arrow_array
                                                             .buffers[2]);
                        const uint32_t*
                            offsets = (const uint32_t*)(arrow_array.buffers[1]);

                        for (int32_t i = 0; i < arrow_array.length; i++) {
                            auto value = std::string{
                                data + offsets[i], offsets[i + 1] - offsets[i]};
                            reader.set_dim_point(dim, value);
                        }

                    } else if (
                        !strcmp(arrow_schema.format, "tss:") ||
                        !strcmp(arrow_schema.format, "tsm:") ||
                        !strcmp(arrow_schema.format, "tsu:") ||
                        !strcmp(arrow_schema.format, "tsn:")) {
                        tcb::span<int64_t> data{
                            (int64_t*)arrow_array.buffers[data_index],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);

                        // TODO:
                        // (pa.bool_(),) * 2,

                    } else if (
                        !strcmp(arrow_schema.format, "U") ||
                        !strcmp(arrow_schema.format, "Z")) {
                        // TODO: partitioning is not supported for string/bytes
                        // dims
                        const char* data = (const char*)(arrow_array
                                                             .buffers[2]);
                        const uint64_t*
                            offsets = (const uint64_t*)(arrow_array.buffers[1]);

                        for (int64_t i = 0; i < arrow_array.length; i++) {
                            auto value = std::string{
                                data + offsets[i], offsets[i + 1] - offsets[i]};
                            reader.set_dim_point(dim, value);
                        }

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
                    return to_table(*buffers);
                }

                // No data was read, the query is complete, return nullopt
                return std::nullopt;
            })

        .def("nnz", &SOMAArray::nnz, py::call_guard<py::gil_scoped_release>())

        .def_property_readonly("shape", &SOMAArray::shape);
}
}  // namespace tiledbsoma
