/**
 * @file   libtiledbsoma.cc
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

#ifdef BUILD_COMMIT_HASH
#define VERSION BUILD_COMMIT_HASH
#else
#define VERSION "dev"
#endif

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

std::string version() {
    int major, minor, patch;
    tiledb_version(&major, &minor, &patch);
    return fmt::format(
        "libtiledbsoma={}\nlibtiledb={}.{}.{}", VERSION, major, minor, patch);
}

/**
 * @brief pybind11 bindings
 *
 */
PYBIND11_MODULE(libtiledbsoma, m) {
    tiledbpy::init_query_condition(m);

    m.doc() = "SOMA acceleration library";

    m.def("version", []() { return version(); });

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
        "stats_enable",
        []() { tiledb::Stats::enable(); },
        "Enable TileDB internal statistics.");
    m.def(
        "stats_disable",
        []() { tiledb::Stats::disable(); },
        "Disable TileDB internal statistics.");
    m.def(
        "stats_reset",
        []() { tiledb::Stats::reset(); },
        "Reset all TileDB internal statistics to 0.");
    m.def(
        "stats_dump",
        []() {
            std::string stats;
            tiledb::Stats::dump(&stats);
            std::cout << version() << "\n" << stats;
        },
        "Print TileDB internal statistics.");

    py::class_<SOMAReader>(m, "SOMAReader")
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

                    auto reader = SOMAReader::open(
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
            [](SOMAReader& reader,
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

                // Reset state of the existing SOMAReader object
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

        // Binding overloaded methods to templated member functions requires
        // more effort, see:
        // https://pybind11.readthedocs.io/en/stable/classes.html#overloaded-methods
        .def(
            "set_dim_points",
            static_cast<void (SOMAReader::*)(
                const std::string&, const std::vector<int64_t>&)>(
                &SOMAReader::set_dim_points))

        .def(
            "set_dim_points",
            static_cast<void (SOMAReader::*)(
                const std::string&, const std::vector<std::string>&)>(
                &SOMAReader::set_dim_points))

        // Binding to set slices using PyArrow::ChunkedArray
        .def(
            "set_dim_points",
            [](SOMAReader& reader,
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
                    array.attr("_export_to_c")(
                        arrow_array_ptr, arrow_schema_ptr);

                    // 'l' = arrow data type int64
                    if (!strcmp(arrow_schema.format, "l")) {
                        tcb::span<int64_t> data{
                            (int64_t*)arrow_array.buffers[1],
                            (uint64_t)arrow_array.length};
                        reader.set_dim_points(
                            dim, data, partition_index, partition_count);
                    } else if (!strcmp(arrow_schema.format, "U")) {
                        // TODO: partitioning is not supported for string dims
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
                            "[libtiledbsoma] set_dim_points: type={} not "
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

        .def(
            "set_dim_ranges",
            static_cast<void (SOMAReader::*)(
                const std::string&,
                const std::vector<std::pair<int64_t, int64_t>>&)>(
                &SOMAReader::set_dim_ranges))

        .def(
            "set_dim_ranges",
            static_cast<void (SOMAReader::*)(
                const std::string&,
                const std::vector<std::pair<std::string, std::string>>&)>(
                &SOMAReader::set_dim_ranges))

        .def(
            "submit",
            &SOMAReader::submit,
            py::call_guard<py::gil_scoped_release>())

        .def("results_complete", &SOMAReader::results_complete)

        .def(
            "read_next",
            [](SOMAReader& reader) -> std::optional<py::object> {
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

        .def("nnz", &SOMAReader::nnz, py::call_guard<py::gil_scoped_release>());
}
}  // namespace tiledbsoma
