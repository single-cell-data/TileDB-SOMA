/**
 * @file   soma_column.cc
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

void load_soma_column(py::module& m) {
    py::class_<SOMAColumn, std::shared_ptr<SOMAColumn>>(m.attr("SOMAColumn"))
        .def(
            "select_columns",
            [](std::shared_ptr<SOMAColumn>& column, ManagedQuery& query) {
                column->select_columns(query);
            })
        .def(
            "set_dim_points_string_or_bytes",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::string>& points) {
                column->set_dim_points<std::string>(mq, points);
            })
        .def(
            "set_dim_points_double",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<double_t>& points) {
                column->set_dim_points<double_t>(mq, points);
            })
        .def(
            "set_dim_points_float",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<float_t>& points) {
                column->set_dim_points<float_t>(mq, points);
            })
        .def(
            "set_dim_points_int64",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<int64_t>& points) {
                column->set_dim_points<int64_t>(mq, points);
            })
        .def(
            "set_dim_points_int32",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<int32_t>& points) {
                column->set_dim_points<int32_t>(mq, points);
            })
        .def(
            "set_dim_points_int16",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<int16_t>& points) {
                column->set_dim_points<int16_t>(mq, points);
            })
        .def(
            "set_dim_points_int8",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<int8_t>& points) {
                column->set_dim_points<int8_t>(mq, points);
            })
        .def(
            "set_dim_points_uint64",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<uint64_t>& points) {
                column->set_dim_points<uint64_t>(mq, points);
            })
        .def(
            "set_dim_points_uint32",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<uint32_t>& points) {
                column->set_dim_points<uint32_t>(mq, points);
            })
        .def(
            "set_dim_points_uint16",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<uint16_t>& points) {
                column->set_dim_points<uint16_t>(mq, points);
            })
        .def(
            "set_dim_points_uint8",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<uint8_t>& points) {
                column->set_dim_points<uint8_t>(mq, points);
            })
        .def(
            "set_dim_points_double_array",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::vector<double_t>>& points) {
                column->set_dim_points<std::vector<double_t>>(mq, points);
            })
        .def(
            "set_dim_ranges_string_or_bytes",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<std::string, std::string>>& ranges) {
                column->set_dim_ranges<std::string>(mq, ranges);
            })
        .def(
            "set_dim_ranges_double",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<double_t, double_t>>& ranges) {
                column->set_dim_ranges<double_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_float",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<float_t, float_t>>& ranges) {
                column->set_dim_ranges<float_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_int64",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<int64_t, int64_t>>& ranges) {
                column->set_dim_ranges<int64_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_int32",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<int32_t, int32_t>>& ranges) {
                column->set_dim_ranges<int32_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_int16",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<int16_t, int16_t>>& ranges) {
                column->set_dim_ranges<int16_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_int8",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<int8_t, int8_t>>& ranges) {
                column->set_dim_ranges<int8_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_uint64",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<uint64_t, uint64_t>>& ranges) {
                column->set_dim_ranges<uint64_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_uint32",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<uint32_t, uint32_t>>& ranges) {
                column->set_dim_ranges<uint32_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_uint16",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<uint16_t, uint16_t>>& ranges) {
                column->set_dim_ranges<uint16_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_uint8",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<std::pair<uint8_t, uint8_t>>& ranges) {
                column->set_dim_ranges<uint8_t>(mq, ranges);
            })
        .def(
            "set_dim_ranges_double_array",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
               const std::vector<
                   std::pair<std::vector<double_t>, std::vector<double_t>>>&
                   ranges) {
                column->set_dim_ranges<std::vector<double_t>>(mq, ranges);
            })
        .def(
            "set_dim_points_arrow",
            [](std::shared_ptr<SOMAColumn>& column,
               ManagedQuery& mq,
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
                            column->set_dim_points<int64_t>(
                                mq, coords.cast<std::vector<int64_t>>());
                        } else if (!strcmp(arrow_schema.format, "i")) {
                            column->set_dim_points<int32_t>(
                                mq, coords.cast<std::vector<int32_t>>());
                        } else if (!strcmp(arrow_schema.format, "s")) {
                            column->set_dim_points<int16_t>(
                                mq, coords.cast<std::vector<int16_t>>());
                        } else if (!strcmp(arrow_schema.format, "c")) {
                            column->set_dim_points<int8_t>(
                                mq, coords.cast<std::vector<int8_t>>());
                        } else if (!strcmp(arrow_schema.format, "L")) {
                            column->set_dim_points<uint64_t>(
                                mq, coords.cast<std::vector<uint64_t>>());
                        } else if (!strcmp(arrow_schema.format, "I")) {
                            column->set_dim_points<uint32_t>(
                                mq, coords.cast<std::vector<uint32_t>>());
                        } else if (!strcmp(arrow_schema.format, "S")) {
                            column->set_dim_points<uint16_t>(
                                mq, coords.cast<std::vector<uint16_t>>());
                        } else if (!strcmp(arrow_schema.format, "C")) {
                            column->set_dim_points<uint8_t>(
                                mq, coords.cast<std::vector<uint8_t>>());
                        } else if (!strcmp(arrow_schema.format, "f")) {
                            column->set_dim_points<float_t>(
                                mq, coords.cast<std::vector<float_t>>());
                        } else if (!strcmp(arrow_schema.format, "g")) {
                            column->set_dim_points<double_t>(
                                mq, coords.cast<std::vector<double_t>>());
                        } else if (
                            !strcmp(arrow_schema.format, "u") ||
                            !strcmp(arrow_schema.format, "U") ||
                            !strcmp(arrow_schema.format, "z") ||
                            !strcmp(arrow_schema.format, "Z")) {
                            column->set_dim_points<std::string>(
                                mq, coords.cast<std::vector<std::string>>());
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
                            column->set_dim_points<int64_t>(
                                mq, coords.cast<std::vector<int64_t>>());
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
            "mq"_a,
            "py_arrow_array"_a,
            "partition_index"_a = 0,
            "partition_count"_a = 1);
}
}  // namespace libtiledbsomacpp
