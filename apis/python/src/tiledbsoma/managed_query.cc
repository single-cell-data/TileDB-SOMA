/**
 * @file   managed_query.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the ManagedQuery bindings.
 */

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

void load_managed_query(py::module& m) {
    py::class_<ManagedQuery>(m, "ManagedQuery")
        .def(
            py::init([](SOMAArray array,
                        std::shared_ptr<SOMAContext> ctx,
                        std::string_view name) {
                return ManagedQuery(array, ctx->tiledb_ctx(), name);
            }),
            py::arg("array"),
            py::arg("ctx"),
            py::arg("name") = "unnamed")

        .def("is_empty_query", &ManagedQuery::is_empty_query)
        .def("is_complete", &ManagedQuery::is_complete)

        .def("set_layout", &ManagedQuery::set_layout)
        .def(
            "set_condition",
            [](ManagedQuery& mq,
               py::object py_query_condition,
               py::object py_schema) {
                auto column_names = mq.column_names();
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
                mq.select_columns(column_names, false, true);

                // Release python GIL after we're done accessing python
                // objects
                py::gil_scoped_release release;
                // Set query condition if present
                if (qc) {
                    mq.set_condition(*qc);
                }
            },
            "py_query_condition"_a,
            "py_schema"_a)
        .def(
            "select_columns",
            &ManagedQuery::select_columns,
            "names"_a,
            "if_not_empty"_a = false,
            "replace"_a = false)

        .def(
            "next",
            [](ManagedQuery& mq) -> std::optional<py::object> {
                // Release python GIL before reading data
                py::gil_scoped_release release;
                std::optional<std::shared_ptr<ArrayBuffers>> tbl;
                try {
                    tbl = mq.read_next();
                    // Acquire python GIL before accessing python objects
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
                py::gil_scoped_acquire acquire;

                if (!tbl) {
                    throw py::stop_iteration();
                }

                return to_table(tbl);
            })

        .def(
            "set_array_data",
            [](ManagedQuery& mq, py::handle py_batch) {
                ArrowSchema arrow_schema;
                ArrowArray arrow_array;
                uintptr_t arrow_schema_ptr = (uintptr_t)(&arrow_schema);
                uintptr_t arrow_array_ptr = (uintptr_t)(&arrow_array);
                py_batch.attr("_export_to_c")(
                    arrow_array_ptr, arrow_schema_ptr);

                py::gil_scoped_release release;
                try {
                    mq.set_array_data(
                        std::make_unique<ArrowSchema>(arrow_schema),
                        std::make_unique<ArrowArray>(arrow_array));
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
                py::gil_scoped_acquire acquire;

                arrow_schema.release(&arrow_schema);
                arrow_array.release(&arrow_array);
            })
        .def(
            "set_column_data",
            [](ManagedQuery& mq, std::string name, py::array data) {
                py::buffer_info data_info = data.request();

                py::gil_scoped_release release;
                try {
                    mq.setup_write_column(
                        name,
                        data.size(),
                        (const void*)data_info.ptr,
                        static_cast<uint64_t*>(nullptr),
                        std::nullopt);
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
                py::gil_scoped_acquire acquire;
            })
        .def(
            "submit_write",
            [](ManagedQuery& mq, bool sort_coords) {
                try {
                    mq.submit_write(sort_coords);
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            },
            "sort_coords"_a = false,
            py::call_guard<py::gil_scoped_release>())

        .def("reset", &ManagedQuery::reset)
        .def("close", &ManagedQuery::close)

        .def_property_readonly("result_order", &ManagedQuery::result_order)
        .def_property_readonly("column_names", &ManagedQuery::column_names)

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
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::string>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_double",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<double_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_float",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<float_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_int64",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<int64_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_int32",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<int32_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_int16",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<int16_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_int8",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<int8_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_uint64",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<uint64_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_uint32",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<uint32_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_uint16",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<uint16_t>& points) {
                try {
                    mq.select_points(dim, points);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_points_uint8",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<uint8_t>& points) {
                try {
                    mq.select_points(dim, points);
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
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<std::string, std::string>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_double",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<double, double>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_float",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<float, float>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_int64",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<int64_t, int64_t>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_int32",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<int32_t, int32_t>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_int16",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<int16_t, int16_t>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_int8",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<int8_t, int8_t>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_uint64",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<uint64_t, uint64_t>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_uint32",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<uint32_t, uint32_t>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_uint16",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<uint16_t, uint16_t>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        .def(
            "set_dim_ranges_uint8",
            [](ManagedQuery& mq,
               const std::string& dim,
               const std::vector<std::pair<uint8_t, uint8_t>>& ranges) {
                try {
                    mq.select_ranges(dim, ranges);
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            })

        // After this are short functions expected to be invoked when the coords
        // are Python list/tuple, or NumPy arrays.  Arrow arrays are in this
        // long if-else-if function.
        .def(
            "set_dim_points_arrow",
            [](ManagedQuery& mq,
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
                            mq.select_points(
                                dim, coords.cast<std::vector<int64_t>>());
                        } else if (!strcmp(arrow_schema.format, "i")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<int32_t>>());
                        } else if (!strcmp(arrow_schema.format, "s")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<int16_t>>());
                        } else if (!strcmp(arrow_schema.format, "c")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<int8_t>>());
                        } else if (!strcmp(arrow_schema.format, "L")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<uint64_t>>());
                        } else if (!strcmp(arrow_schema.format, "I")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<uint32_t>>());
                        } else if (!strcmp(arrow_schema.format, "S")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<uint16_t>>());
                        } else if (!strcmp(arrow_schema.format, "C")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<uint8_t>>());
                        } else if (!strcmp(arrow_schema.format, "f")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<float>>());
                        } else if (!strcmp(arrow_schema.format, "g")) {
                            mq.select_points(
                                dim, coords.cast<std::vector<double>>());
                        } else if (
                            !strcmp(arrow_schema.format, "u") ||
                            !strcmp(arrow_schema.format, "z")) {
                            mq.select_points(
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
                            mq.select_points(
                                dim, coords.cast<std::vector<int64_t>>());
                        } else if (
                            !strcmp(arrow_schema.format, "U") ||
                            !strcmp(arrow_schema.format, "Z")) {
                            mq.select_points(
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
            "partition_count"_a = 1);
}
}  // namespace libtiledbsomacpp
