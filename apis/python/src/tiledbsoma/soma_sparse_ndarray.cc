/**
 * @file   soma_sparse_ndarray.cc
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
 * This file defines the SOMASparseNDArray bindings.
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

void load_soma_sparse_ndarray(py::module& m) {
    py::class_<SOMASparseNDArray, SOMAArray, SOMAObject>(m, "SOMASparseNDArray")

        .def_static(
            "create",
            [](std::string_view uri,
               std::string format,
               py::object index_column_info,
               std::shared_ptr<SOMAContext> context,
               PlatformConfig platform_config,
               std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                ArrowSchema index_column_schema;
                ArrowArray index_column_array;
                uintptr_t
                    index_column_schema_ptr = (uintptr_t)(&index_column_schema);
                uintptr_t
                    index_column_array_ptr = (uintptr_t)(&index_column_array);
                index_column_info.attr("_export_to_c")(
                    index_column_array_ptr, index_column_schema_ptr);

                try {
                    SOMASparseNDArray::create(
                        uri,
                        format,
                        ArrowTable(
                            std::make_unique<ArrowArray>(index_column_array),
                            std::make_unique<ArrowSchema>(index_column_schema)),
                        context,
                        platform_config,
                        timestamp);
                } catch (const std::out_of_range& e) {
                    throw py::type_error(e.what());
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
                index_column_array.release(&index_column_array);
                index_column_schema.release(&index_column_schema);
            },
            "uri"_a,
            py::kw_only(),
            "format"_a,
            "index_column_info"_a,
            "ctx"_a,
            "platform_config"_a,
            "timestamp"_a = py::none())

        .def_static(
            "open",
            py::overload_cast<
                std::string_view,
                OpenMode,
                std::shared_ptr<SOMAContext>,
                std::vector<std::string>,
                ResultOrder,
                std::optional<std::pair<uint64_t, uint64_t>>>(
                &SOMASparseNDArray::open),
            "uri"_a,
            "mode"_a,
            "context"_a,
            py::kw_only(),
            "column_names"_a = py::tuple(),
            "result_order"_a = ResultOrder::automatic,
            "timestamp"_a = py::none())

        .def_static("exists", &SOMASparseNDArray::exists)

        .def_property_readonly("shape", &SOMASparseNDArray::shape)
        .def_property_readonly("maxshape", &SOMASparseNDArray::maxshape)
        .def_property_readonly(
            "tiledbsoma_has_upgraded_shape", &SOMAArray::has_current_domain)

        .def(
            "resize",
            [](SOMAArray& array, const std::vector<int64_t>& newshape) {
                try {
                    array.resize(newshape, "resize");
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a)

        .def(
            "tiledbsoma_upgrade_shape",
            [](SOMAArray& array, const std::vector<int64_t>& newshape) {
                try {
                    array.upgrade_shape(newshape, "tiledbsoma_upgrade_shape");
                } catch (const std::exception& e) {
                    throw TileDBSOMAError(e.what());
                }
            },
            "newshape"_a);
}
}  // namespace libtiledbsomacpp
