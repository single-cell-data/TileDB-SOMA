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

void write(SOMAArray& array, py::handle py_batch, bool sort_coords = true) {
    ArrowSchema arrow_schema;
    ArrowArray arrow_array;
    uintptr_t arrow_schema_ptr = (uintptr_t)(&arrow_schema);
    uintptr_t arrow_array_ptr = (uintptr_t)(&arrow_array);
    py_batch.attr("_export_to_c")(arrow_array_ptr, arrow_schema_ptr);

    try {
        array.set_array_data(
            std::make_unique<ArrowSchema>(arrow_schema),
            std::make_unique<ArrowArray>(arrow_array));
    } catch (const std::exception& e) {
        TPY_ERROR_LOC(e.what());
    }

    try {
        array.write(sort_coords);
    } catch (const std::exception& e) {
        TPY_ERROR_LOC(e.what());
    }
    arrow_schema.release(&arrow_schema);
    arrow_array.release(&arrow_array);
}

void write_coords(
    SOMAArray& array,
    std::vector<py::array> coords,
    py::array data,
    bool sort_coords = true) {
    for (uint64_t i = 0; i < coords.size(); ++i) {
        py::buffer_info coords_info = coords[i].request();
        array.set_column_data(
            "soma_dim_" + std::to_string(i),
            coords[i].size(),
            (const void*)coords_info.ptr);
    }

    py::buffer_info data_info = data.request();
    array.set_column_data("soma_data", data.size(), (const void*)data_info.ptr);

    try {
        array.write(sort_coords);
    } catch (const std::exception& e) {
        TPY_ERROR_LOC(e.what());
    }
}

// This is shared code for non-empty domain, domain, and maxdomain.  Returns a
// Python list of Arrow arrays. Python code can convert these further.
py::list domainish_to_list(ArrowArray* arrow_array, ArrowSchema* arrow_schema) {
    auto pa = py::module::import("pyarrow");
    auto pa_array_import = pa.attr("Array").attr("_import_from_c");

    py::list array_list;
    for (int i = 0; i < arrow_array->n_children; i++) {
        // Note that this runs the release callbacks within the ArrowArray and
        // the ArrowSchema, freeing memory.
        auto array = pa_array_import(
            py::capsule(arrow_array->children[i]),
            py::capsule(arrow_schema->children[i]));
        array_list.append(array);

        // Already released: ensure there is no attempt at second free.
        arrow_array->children[i] = nullptr;
        arrow_schema->children[i] = nullptr;
    }
    // Already released: ensure there is no attempt at second free.
    arrow_array->n_children = 0;
    arrow_array->children = nullptr;
    arrow_schema->n_children = 0;
    arrow_schema->children = nullptr;

    return array_list;
}

void load_soma_array(py::module& m) {
    py::class_<SOMAArray, SOMAObject>(m, "SOMAArray")
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
                        std::make_shared<SOMAContext>(platform_config),
                        name,
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

        .def("__enter__", [](SOMAArray& array) { return array; })
        .def(
            "__exit__",
            [](SOMAArray& array,
               py::object exc_type,
               py::object exc_value,
               py::object traceback) { array.close(); })

        .def(
            "set_condition",
            [](SOMAArray& array,
               py::object py_query_condition,
               py::object py_schema) {
                auto column_names = array.column_names();
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
                array.reset(column_names);

                // Release python GIL after we're done accessing python
                // objects
                py::gil_scoped_release release;
                // Set query condition if present
                if (qc) {
                    array.set_condition(*qc);
                }
            },
            "py_query_condition"_a,
            "py_schema"_a)

        .def(
            "reset",
            [](SOMAArray& array,
               std::optional<std::vector<std::string>> column_names_in,
               std::string_view batch_size,
               ResultOrder result_order) {
                // Handle optional args
                std::vector<std::string> column_names;
                if (column_names_in) {
                    column_names = *column_names_in;
                }

                // Reset state of the existing SOMAArray object
                array.reset(column_names, batch_size, result_order);
            },
            py::kw_only(),
            "column_names"_a = py::none(),
            "batch_size"_a = "auto",
            "result_order"_a = ResultOrder::automatic)

        .def("reopen", &SOMAArray::reopen)
        .def("close", &SOMAArray::close)
        .def_property_readonly(
            "closed",
            [](SOMAArray& array) -> bool { return not array.is_open(); })
        .def_property_readonly(
            "mode",
            [](SOMAArray& array) {
                return array.mode() == OpenMode::read ? "r" : "w";
            })
        .def_property_readonly(
            "schema",
            [](SOMAArray& array) -> py::object {
                auto pa = py::module::import("pyarrow");
                auto pa_schema_import = pa.attr("Schema").attr(
                    "_import_from_c");
                return pa_schema_import(
                    py::capsule(array.arrow_schema().get()));
            })
        .def("context", &SOMAArray::ctx)

        // After this are short functions expected to be invoked when the coords
        // are Python list/tuple, or NumPy arrays.  Arrow arrays are in this
        // long if-else-if function.
        .def(
            "set_dim_points_arrow",
            [](SOMAArray& array,
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

                for (const pybind11::handle array_handle : array_chunks) {
                    ArrowSchema arrow_schema;
                    ArrowArray arrow_array;
                    uintptr_t arrow_schema_ptr = (uintptr_t)(&arrow_schema);
                    uintptr_t arrow_array_ptr = (uintptr_t)(&arrow_array);

                    // Call handle._export_to_c to get arrow array and schema
                    //
                    // If ever a NumPy array gets in here, there will be an
                    // exception like "AttributeError: 'numpy.ndarray' object
                    // has no attribute '_export_to_c'".
                    array_handle.attr("_export_to_c")(
                        arrow_array_ptr, arrow_schema_ptr);

                    auto coords = array_handle.attr("tolist")();

                    try {
                        if (!strcmp(arrow_schema.format, "l")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<int64_t>>());
                        } else if (!strcmp(arrow_schema.format, "i")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<int32_t>>());
                        } else if (!strcmp(arrow_schema.format, "s")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<int16_t>>());
                        } else if (!strcmp(arrow_schema.format, "c")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<int8_t>>());
                        } else if (!strcmp(arrow_schema.format, "L")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<uint64_t>>());
                        } else if (!strcmp(arrow_schema.format, "I")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<uint32_t>>());
                        } else if (!strcmp(arrow_schema.format, "S")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<uint16_t>>());
                        } else if (!strcmp(arrow_schema.format, "C")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<uint8_t>>());
                        } else if (!strcmp(arrow_schema.format, "f")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<float>>());
                        } else if (!strcmp(arrow_schema.format, "g")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<double>>());
                        } else if (
                            !strcmp(arrow_schema.format, "u") ||
                            !strcmp(arrow_schema.format, "z")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<std::string>>());
                        } else if (
                            !strcmp(arrow_schema.format, "tss:") ||
                            !strcmp(arrow_schema.format, "tsm:") ||
                            !strcmp(arrow_schema.format, "tsu:") ||
                            !strcmp(arrow_schema.format, "tsn:")) {
                            // convert the Arrow Array to int64
                            auto pa = py::module::import("pyarrow");
                            coords = array_handle
                                         .attr("cast")(pa.attr("int64")())
                                         .attr("tolist")();
                            array.set_dim_points(
                                dim, coords.cast<std::vector<int64_t>>());
                        } else if (
                            !strcmp(arrow_schema.format, "U") ||
                            !strcmp(arrow_schema.format, "Z")) {
                            array.set_dim_points(
                                dim, coords.cast<std::vector<std::string>>());
                        } else {
                            TPY_ERROR_LOC(
                                "[pytiledbsoma] set_dim_points: type={} not "
                                "supported" +
                                std::string(arrow_schema.format));
                        }
                    } catch (const std::exception& e) {
                        throw TileDBSOMAError(e.what());
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
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::string>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_double",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<double>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_float",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<float>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_int64",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<int64_t>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_int32",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<int32_t>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_int16",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<int16_t>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_int8",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<int8_t>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_uint64",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<uint64_t>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_uint32",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<uint32_t>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_uint16",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<uint16_t>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_uint8",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<uint8_t>& points) {
                try {
                    array.set_dim_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

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
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<std::string, std::string>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_double",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<double, double>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_float",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<float, float>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_int64",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<int64_t, int64_t>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_int32",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<int32_t, int32_t>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_int16",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<int16_t, int16_t>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_int8",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<int8_t, int8_t>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_uint64",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<uint64_t, uint64_t>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_uint32",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<uint32_t, uint32_t>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_uint16",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<uint16_t, uint16_t>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_uint8",
            [](SOMAArray& array,
               const std::string& dim,
               const std::vector<std::pair<uint8_t, uint8_t>>& ranges) {
                try {
                    array.set_dim_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def("results_complete", &SOMAArray::results_complete)

        .def(
            "read_next",
            [](SOMAArray& array) -> std::optional<py::object> {
                // Release python GIL before reading data
                py::gil_scoped_release release;

                // Try to read more data
                try {
                    auto buffers = array.read_next();

                    // If more data was read, convert it to an arrow table and
                    // return
                    if (buffers.has_value()) {
                        // Acquire python GIL before accessing python objects
                        py::gil_scoped_acquire acquire;
                        return to_table(*buffers);
                    }
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }

                // No data was read, the query is complete, return nullopt
                return std::nullopt;
            })

        .def("write", write)

        .def("write_coords", write_coords)

        .def("nnz", &SOMAArray::nnz, py::call_guard<py::gil_scoped_release>())

        .def_property_readonly("uri", &SOMAArray::uri)

        .def_property_readonly("column_names", &SOMAArray::column_names)

        .def_property_readonly("result_order", &SOMAArray::result_order)

        .def_property_readonly(
            "timestamp",
            [](SOMAArray& array) -> py::object {
                if (!array.timestamp().has_value())
                    return py::none();
                return py::cast(array.timestamp()->second);
            })

        .def("domainish_to_list", domainish_to_list)

        .def(
            "non_empty_domain",
            [](SOMAArray& array) {
                ArrowTable arrow_table = array.get_non_empty_domain();
                return domainish_to_list(
                    arrow_table.first.get(), arrow_table.second.get());
            })

        .def(
            "domain",
            [](SOMAArray& array) {
                auto pa = py::module::import("pyarrow");
                ArrowTable arrow_table = array.get_soma_domain();
                return domainish_to_list(
                    arrow_table.first.get(), arrow_table.second.get());
            })

        .def(
            "maxdomain",
            [](SOMAArray& array) {
                auto pa = py::module::import("pyarrow");
                ArrowTable arrow_table = array.get_soma_maxdomain();
                return domainish_to_list(
                    arrow_table.first.get(), arrow_table.second.get());
            })

        .def(
            "non_empty_domain_slot",
            [](SOMAArray& array, std::string name, py::dtype dtype) {
                switch (np_to_tdb_dtype(dtype)) {
                    case TILEDB_UINT64:
                        return py::cast(
                            array.non_empty_domain_slot<uint64_t>(name));
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
                        return py::cast(
                            array.non_empty_domain_slot<int64_t>(name));
                    case TILEDB_UINT32:
                        return py::cast(
                            array.non_empty_domain_slot<uint32_t>(name));
                    case TILEDB_INT32:
                        return py::cast(
                            array.non_empty_domain_slot<int32_t>(name));
                    case TILEDB_UINT16:
                        return py::cast(
                            array.non_empty_domain_slot<uint16_t>(name));
                    case TILEDB_INT16:
                        return py::cast(
                            array.non_empty_domain_slot<int16_t>(name));
                    case TILEDB_UINT8:
                        return py::cast(
                            array.non_empty_domain_slot<uint8_t>(name));
                    case TILEDB_INT8:
                        return py::cast(
                            array.non_empty_domain_slot<int8_t>(name));
                    case TILEDB_FLOAT64:
                        return py::cast(
                            array.non_empty_domain_slot<double>(name));
                    case TILEDB_FLOAT32:
                        return py::cast(
                            array.non_empty_domain_slot<float>(name));
                    case TILEDB_STRING_UTF8:
                    case TILEDB_STRING_ASCII:
                        return py::cast(array.non_empty_domain_slot_var(name));
                    default:
                        throw TileDBSOMAError(
                            "Unsupported dtype for nonempty domain.");
                }
            })

        .def(
            "soma_domain_slot",
            [](SOMAArray& array, std::string name, py::dtype dtype) {
                switch (np_to_tdb_dtype(dtype)) {
                    case TILEDB_UINT64:
                        return py::cast(array.soma_domain_slot<uint64_t>(name));
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
                        return py::cast(array.soma_domain_slot<int64_t>(name));
                    case TILEDB_UINT32:
                        return py::cast(array.soma_domain_slot<uint32_t>(name));
                    case TILEDB_INT32:
                        return py::cast(array.soma_domain_slot<int32_t>(name));
                    case TILEDB_UINT16:
                        return py::cast(array.soma_domain_slot<uint16_t>(name));
                    case TILEDB_INT16:
                        return py::cast(array.soma_domain_slot<int16_t>(name));
                    case TILEDB_UINT8:
                        return py::cast(array.soma_domain_slot<uint8_t>(name));
                    case TILEDB_INT8:
                        return py::cast(array.soma_domain_slot<int8_t>(name));
                    case TILEDB_FLOAT64:
                        return py::cast(array.soma_domain_slot<double>(name));
                    case TILEDB_FLOAT32:
                        return py::cast(array.soma_domain_slot<float>(name));
                    case TILEDB_STRING_UTF8:
                    case TILEDB_STRING_ASCII: {
                        std::pair<std::string, std::string> str_domain;
                        return py::cast(std::make_pair("", ""));
                    }
                    default:
                        throw TileDBSOMAError(
                            "Unsupported dtype for Dimension's domain");
                }
            })

        .def(
            "soma_maxdomain_slot",
            [](SOMAArray& array, std::string name, py::dtype dtype) {
                switch (np_to_tdb_dtype(dtype)) {
                    case TILEDB_UINT64:
                        return py::cast(
                            array.soma_maxdomain_slot<uint64_t>(name));
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
                        return py::cast(
                            array.soma_maxdomain_slot<int64_t>(name));
                    case TILEDB_UINT32:
                        return py::cast(
                            array.soma_maxdomain_slot<uint32_t>(name));
                    case TILEDB_INT32:
                        return py::cast(
                            array.soma_maxdomain_slot<int32_t>(name));
                    case TILEDB_UINT16:
                        return py::cast(
                            array.soma_maxdomain_slot<uint16_t>(name));
                    case TILEDB_INT16:
                        return py::cast(
                            array.soma_maxdomain_slot<int16_t>(name));
                    case TILEDB_UINT8:
                        return py::cast(
                            array.soma_maxdomain_slot<uint8_t>(name));
                    case TILEDB_INT8:
                        return py::cast(
                            array.soma_maxdomain_slot<int8_t>(name));
                    case TILEDB_FLOAT64:
                        return py::cast(
                            array.soma_maxdomain_slot<double>(name));
                    case TILEDB_FLOAT32:
                        return py::cast(array.soma_maxdomain_slot<float>(name));
                    case TILEDB_STRING_UTF8:
                    case TILEDB_STRING_ASCII: {
                        std::pair<std::string, std::string> str_domain;
                        return py::cast(std::make_pair("", ""));
                    }
                    default:
                        throw TileDBSOMAError(
                            "Unsupported dtype for Dimension's domain");
                }
            })

        .def_property_readonly("dimension_names", &SOMAArray::dimension_names)

        .def("consolidate_and_vacuum", &SOMAArray::consolidate_and_vacuum)

        .def_property_readonly(
            "meta",
            [](SOMAArray& array) -> py::dict {
                return meta(array.get_metadata());
            })

        .def("set_metadata", set_metadata)

        .def("delete_metadata", &SOMAArray::delete_metadata)

        .def(
            "get_metadata",
            py::overload_cast<const std::string&>(&SOMAArray::get_metadata))

        .def("has_metadata", &SOMAArray::has_metadata)

        .def("metadata_num", &SOMAArray::metadata_num);
}
}  // namespace libtiledbsomacpp
