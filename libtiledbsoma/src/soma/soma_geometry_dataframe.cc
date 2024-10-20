/**
 * @file   soma_geometry_dataframe.cc
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
 *   This file defines the SOMAGeometryDataFrame class.
 */

#include "soma_geometry_dataframe.h"
#include "../utils/util.h"

#include <regex>

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAGeometryDataFrame::create(
    std::string_view uri,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    ArrowTable spatial_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    std::vector<std::string> spatial_axes;
    auto tiledb_schema = ArrowAdapter::tiledb_schema_from_arrow_schema(
        ctx->tiledb_ctx(),
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ArrowTable(
            std::move(spatial_columns.first),
            std::move(spatial_columns.second)),
        "SOMAGeometryDataFrame",
        true,
        platform_config);
    auto array = SOMAArray::create(
        ctx, uri, tiledb_schema, "SOMAGeometryDataFrame", timestamp);
}

std::unique_ptr<SOMAGeometryDataFrame> SOMAGeometryDataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp) {
    return std::make_unique<SOMAGeometryDataFrame>(
        mode, uri, ctx, column_names, result_order, timestamp);
}

bool SOMAGeometryDataFrame::exists(
    std::string_view uri, std::shared_ptr<SOMAContext> ctx) {
    try {
        auto obj = SOMAObject::open(uri, OpenMode::read, ctx);
        return "SOMAGeometryDataFrame" == obj->type();
    } catch (TileDBSOMAError& e) {
        return false;
    }
}

//===================================================================
//= public non-static
//===================================================================

std::unique_ptr<ArrowSchema> SOMAGeometryDataFrame::schema() const {
    return this->arrow_schema();
}

const std::vector<std::string> SOMAGeometryDataFrame::index_column_names()
    const {
    return this->dimension_names();
}

const std::vector<std::string> SOMAGeometryDataFrame::spatial_column_names()
    const {
    std::vector<std::string> names;
    std::unordered_set<std::string> unique_names;
    std::regex rgx("tiledb__internal__(\\S+)__");
    std::smatch matches;
    for (auto dimension : this->dimension_names()) {
        if (std::regex_search(dimension, matches, rgx)) {
            if (unique_names.count(matches[1].str()) == 0) {
                unique_names.insert(matches[1].str());
                names.push_back(matches[1].str());
            }
        }
    }

    return names;
}

uint64_t SOMAGeometryDataFrame::count() {
    return this->nnz();
}

void SOMAGeometryDataFrame::set_array_data(
    std::unique_ptr<ArrowSchema> arrow_schema,
    std::unique_ptr<ArrowArray> arrow_array) {
    std::vector<std::string> spatial_axes = this->spatial_column_names();

    for (auto i = 0; i < arrow_schema->n_children; ++i) {
        if (strcmp(arrow_schema->children[i]->name, "soma_geometry") == 0 &&
            strcmp(arrow_schema->children[i]->format, "+l") == 0) {
            std::vector<ArrowArray*> arrays = util::cast_vertices_to_wkb(
                arrow_array->children[i], spatial_axes);

            // Reconstruct new array and schema
            std::unique_ptr<ArrowSchema>
                transformed_arrow_schema = std::make_unique<ArrowSchema>(
                    ArrowSchema{});
            std::unique_ptr<ArrowArray>
                transformed_arrow_array = std::make_unique<ArrowArray>(
                    ArrowArray{});

            NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
                transformed_arrow_schema.get(),
                ArrowType::NANOARROW_TYPE_STRUCT));
            NANOARROW_THROW_NOT_OK(ArrowSchemaAllocateChildren(
                transformed_arrow_schema.get(),
                arrow_schema->n_children + 2 * spatial_axes.size()));

            NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
                transformed_arrow_array.get(),
                ArrowType::NANOARROW_TYPE_STRUCT));
            NANOARROW_THROW_NOT_OK(ArrowArrayAllocateChildren(
                transformed_arrow_array.get(),
                arrow_array->n_children + 2 * spatial_axes.size()));

            for (int64_t j = 0; j < arrow_schema->n_children; ++j) {
                if (strcmp(arrow_schema->children[j]->name, "soma_geometry") ==
                    0) {
                    NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
                        transformed_arrow_schema->children[j],
                        ArrowType::NANOARROW_TYPE_LARGE_BINARY));
                    NANOARROW_THROW_NOT_OK(ArrowSchemaSetName(
                        transformed_arrow_schema->children[j],
                        "soma_geometry"));

                    ArrowArrayMove(
                        arrays.front(), transformed_arrow_array->children[j]);
                } else {
                    ArrowSchemaMove(
                        arrow_schema->children[j],
                        transformed_arrow_schema->children[j]);
                    ArrowArrayMove(
                        arrow_array->children[j],
                        transformed_arrow_array->children[j]);
                }
            }

            for (size_t j = 0; j < spatial_axes.size(); ++j) {
                NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
                    transformed_arrow_schema
                        ->children[arrow_schema->n_children + 2 * j],
                    ArrowType::NANOARROW_TYPE_DOUBLE));
                NANOARROW_THROW_NOT_OK(ArrowSchemaSetName(
                    transformed_arrow_schema
                        ->children[arrow_schema->n_children + 2 * j],
                    ("tiledb__internal__" + spatial_axes[j] + "__min")
                        .c_str()));
                ArrowArrayMove(
                    arrays[2 * j + 1],
                    transformed_arrow_array
                        ->children[arrow_schema->n_children + 2 * j]);

                NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
                    transformed_arrow_schema
                        ->children[arrow_schema->n_children + 2 * j + 1],
                    ArrowType::NANOARROW_TYPE_DOUBLE));
                NANOARROW_THROW_NOT_OK(ArrowSchemaSetName(
                    transformed_arrow_schema
                        ->children[arrow_schema->n_children + 2 * j + 1],
                    ("tiledb__internal__" + spatial_axes[j] + "__max")
                        .c_str()));
                ArrowArrayMove(
                    arrays[2 * j + 2],
                    transformed_arrow_array
                        ->children[arrow_schema->n_children + 2 * j + 1]);
            }

            SOMAArray::set_array_data(
                std::move(transformed_arrow_schema),
                std::move(transformed_arrow_array));
            return;
        }
    }

    SOMAArray::set_array_data(std::move(arrow_schema), std::move(arrow_array));
}
}  // namespace tiledbsoma
