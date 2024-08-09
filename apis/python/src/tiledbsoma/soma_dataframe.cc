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
               py::object index_column_info,
               std::shared_ptr<SOMAContext> context,
               PlatformConfig platform_config,
               std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
                ArrowSchema schema;
                uintptr_t schema_ptr = (uintptr_t)(&schema);
                py_schema.attr("_export_to_c")(schema_ptr);

                // In the Arrow C API, there is an ArrowSchema for the table.
                // Each of its children is also of type ArrowSchema.
                //
                // Unfortunately, there are two different ways for attributes
                // to be marked nullable:
                //
                // * Flags on the attribute
                //    o Example:
                //      pa.schema(
                //         [pa.field("a", pa.int32())],
                //         [pa.field("b", pa.int32(), nullable=False)],
                //         [pa.field("c", pa.int32(), nullable=True)],
                //     )
                //   o `pa.field` defaults to nullable: here, fields a and c
                //     are both nullable; only b is not.
                //
                // * Metadata for the attribute:
                //     pa.schema(
                //         [pa.field("d", pa.int32())],
                //         metadata={"d": "nullable"},
                //     )
                //
                // Concerns for us:
                //
                // * Arrow fields are nullable by default. So we must make them
                //   non-nullable only when this is explicit in the schema
                //   flags.
                // * If there is no metadata, we simply respect the user's
                //   schema including its per-attribute nullability flags.
                // * If there is set-nullable metadata for an attribute, it
                //   must override the nullability flags for that attribute.

                auto metadata = py_schema.attr("metadata");
                if (py::hasattr(metadata, "get")) {
                    for (int64_t i = 0; i < schema.n_children; ++i) {
                        auto child = schema.children[i];
                        auto val = metadata.attr("get")(
                            py::str(child->name).attr("encode")("utf-8"));

                        if (!val.is(py::none()) &&
                            val.cast<std::string>() == "nullable") {
                            child->flags |= ARROW_FLAG_NULLABLE;
                        }
                    }
                }

                ArrowSchema index_column_schema;
                ArrowArray index_column_array;
                uintptr_t
                    index_column_schema_ptr = (uintptr_t)(&index_column_schema);
                uintptr_t
                    index_column_array_ptr = (uintptr_t)(&index_column_array);
                index_column_info.attr("_export_to_c")(
                    index_column_array_ptr, index_column_schema_ptr);

                try {
                    SOMADataFrame::create(
                        uri,
                        std::make_unique<ArrowSchema>(schema),
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
                schema.release(&schema);
            },
            "uri"_a,
            py::kw_only(),
            "schema"_a,
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
                &SOMADataFrame::open),
            "uri"_a,
            "mode"_a,
            "context"_a,
            py::kw_only(),
            "column_names"_a = py::tuple(),
            "result_order"_a = ResultOrder::automatic,
            "timestamp"_a = py::none())

        .def_static("exists", &SOMADataFrame::exists)
        .def_property_readonly(
            "index_column_names", &SOMADataFrame::index_column_names)
        .def_property_readonly(
            "count",
            &SOMADataFrame::count,
            py::call_guard<py::gil_scoped_release>());
}
}  // namespace libtiledbsomacpp
