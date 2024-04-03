/**
 * @file   soma_dataframe.cc
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
 * This file defines the SOMADataFrame bindings.
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

void load_soma_dataframe(py::module& m) {
    py::class_<SOMADataFrame, SOMAArray, SOMAObject>(m, "SOMADataFrame")

        .def_static(
            "create",
            [](std::string_view uri,
               py::object py_schema,
               std::vector<std::string> index_columns_names,
               py::object py_domains,
               py::object py_extents,
               std::shared_ptr<SOMAContext> context,
               std::optional<PlatformConfig> platform_config,
               std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                ArrowSchema schema;
                uintptr_t schema_ptr = (uintptr_t)(&schema);
                py_schema.attr("_export_to_c")(schema_ptr);

                for (int64_t sch_idx = 0; sch_idx < schema.n_children;
                     ++sch_idx) {
                    auto child = schema.children[sch_idx];
                    auto metadata = py_schema.attr("metadata");
                    if (py::hasattr(metadata, "get")) {
                        auto val = metadata.attr("get")(
                            py::str(child->name).attr("encode")("utf-8"));

                        if (val != py::none() &&
                            val.cast<std::string>() == "nullable") {
                            child->flags &= ARROW_FLAG_NULLABLE;
                        } else {
                            child->flags &= ~ARROW_FLAG_NULLABLE;
                        }
                    }
                }

                ArrowArray domains;
                uintptr_t domains_ptr = (uintptr_t)(&domains);
                py_domains.attr("_export_to_c")(domains_ptr);

                ArrowArray extents;
                uintptr_t extents_ptr = (uintptr_t)(&extents);
                py_extents.attr("_export_to_c")(extents_ptr);

                try {
                    SOMADataFrame::create(
                        uri,
                        std::make_shared<ArrowSchema>(schema),
                        ColumnIndexInfo(
                            index_columns_names,
                            std::make_shared<ArrowArray>(domains),
                            std::make_shared<ArrowArray>(extents)),
                        context,
                        platform_config,
                        timestamp);
                } catch (const std::out_of_range& e) {
                    throw py::type_error(e.what());
                } catch (const std::exception& e) {
                    TPY_ERROR_LOC(e.what());
                }
            })

        .def_static(
            "open",
            py::overload_cast<
                std::string_view,
                OpenMode,
                std::shared_ptr<SOMAContext>,
                std::vector<std::string>,
                ResultOrder,
                std::optional<std::pair<uint64_t, uint64_t>>>(
                &SOMADataFrame::open),
            "uri"_a,
            "mode"_a,
            "context"_a,
            py::kw_only(),
            "column_names"_a = py::none(),
            "result_order"_a = ResultOrder::automatic,
            "timestamp"_a = py::none())

        .def_static("exists", &SOMADataFrame::exists)
        .def_property_readonly(
            "index_column_names", &SOMADataFrame::index_column_names)
        .def_property_readonly("count", &SOMADataFrame::count);
}
}  // namespace libtiledbsomacpp