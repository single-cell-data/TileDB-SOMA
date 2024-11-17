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
#include "../geometry/geometry.h"
#include "../geometry/operators/envelope.h"
#include "../geometry/operators/io/write.h"
#include "../utils/util.h"

#include <regex>
#include <unordered_set>

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAGeometryDataFrame::create(
    std::string_view uri,
    const std::unique_ptr<ArrowSchema>& schema,
    const ArrowTable& index_columns,
    const ArrowTable& spatial_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    std::vector<std::string> spatial_axes;
    auto tiledb_schema = ArrowAdapter::tiledb_schema_from_arrow_schema(
        ctx->tiledb_ctx(),
        schema,
        index_columns,
        "SOMAGeometryDataFrame",
        true,
        platform_config,
        spatial_columns);
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
        /**
         * If `soma_geometry` conforms to specific formats automatically convert
         * to WKB and create additional index columns for spatial axes.
         *
         * If the `soma_geometry` array is a WKB binary users are expected to
         * provide the additional index columns for spatial axes.
         */

        if (strcmp(arrow_schema->children[i]->name, "soma_geometry") == 0 &&
            strcmp(arrow_schema->children[i]->format, "+l") == 0) {
            std::string_view type_metadata;

            if (ArrowMetadataHasKey(
                    arrow_schema->children[i]->metadata,
                    ArrowCharView("geometry_type"))) {
                ArrowStringView out;
                NANOARROW_THROW_NOT_OK(ArrowMetadataGetValue(
                    arrow_schema->children[i]->metadata,
                    ArrowCharView("geometry_type"),
                    &out));

                type_metadata = std::string_view(out.data, out.size_bytes);
            }

            ArrowTable casted_data;
            if (type_metadata == "polygon_ring") {
                auto wkb_data = _cast_polygon_vertex_list_to_wkb(
                    arrow_array->children[i]);
                casted_data = _reconstruct_geometry_data_table(
                    ArrowTable(std::move(arrow_array), std::move(arrow_schema)),
                    wkb_data);
            } else {
                throw std::runtime_error("Unknown geometry type");
            }

            return SOMAArray::set_array_data(
                std::move(casted_data.second), std::move(casted_data.first));
        }
    }

    SOMAArray::set_array_data(std::move(arrow_schema), std::move(arrow_array));
}

//===================================================================
//= private non-static
//===================================================================

std::vector<ArrowTable> SOMAGeometryDataFrame::_cast_polygon_vertex_list_to_wkb(
    ArrowArray* array) {
    // Initialize a vector to hold all the Arrow tables containing the
    // transformed geometry data
    std::vector<std::string> spatial_axes = this->spatial_column_names();
    std::vector<ArrowTable> tables;
    tables.push_back(ArrowTable(
        std::make_unique<ArrowArray>(ArrowArray{}),
        std::make_unique<ArrowSchema>(ArrowSchema{})));

    NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
        tables.front().first.get(), ArrowType::NANOARROW_TYPE_LARGE_BINARY));
    NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
        tables.front().second.get(), ArrowType::NANOARROW_TYPE_LARGE_BINARY));
    NANOARROW_THROW_NOT_OK(
        ArrowSchemaSetName(tables.front().second.get(), "soma_geometry"));

    for (const auto& axis : spatial_axes) {
        // Min spatial axis
        tables.push_back(ArrowTable(
            std::make_unique<ArrowArray>(ArrowArray{}),
            std::make_unique<ArrowSchema>(ArrowSchema{})));
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            tables.back().first.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
            tables.back().second.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaSetName(
            tables.back().second.get(),
            (SOMA_GEOMETRY_DIMENSION_PREFIX + axis + "__min").c_str()));

        // Max spatial axis
        tables.push_back(ArrowTable(
            std::make_unique<ArrowArray>(), std::make_unique<ArrowSchema>()));
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            tables.back().first.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
            tables.back().second.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaSetName(
            tables.back().second.get(),
            (SOMA_GEOMETRY_DIMENSION_PREFIX + axis + "__max").c_str()));
    }

    // Large list of doubles
    const uint32_t* offset = static_cast<const uint32_t*>(array->buffers[1]);
    const double_t* data = static_cast<const double_t*>(
        array->children[0]->buffers[1]);

    size_t wkb_buffer_size = 0;
    std::vector<geometry::GenericGeometry> geometries;

    for (int64_t index = 0; index < array->length; ++index) {
        int64_t stop_index = index < array->length - 1 ?
                                 offset[index + 1] :
                                 array->children[0]->length;

        std::vector<geometry::BasePoint> ring;
        for (int64_t j = offset[index]; j < stop_index; j += 2) {
            ring.push_back(geometry::BasePoint(data[j], data[j + 1]));
        }

        geometries.push_back(
            geometry::GenericGeometry(geometry::Polygon(std::move(ring))));
        wkb_buffer_size += wkb_size(geometries.back());
    }

    NANOARROW_THROW_NOT_OK(
        ArrowArrayReserve(tables.front().first.get(), wkb_buffer_size));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayStartAppending(tables.front().first.get()));
    for (size_t i = 1; i < tables.size(); ++i) {
        NANOARROW_THROW_NOT_OK(
            ArrowArrayReserve(tables[i].first.get(), array->length));
        NANOARROW_THROW_NOT_OK(ArrowArrayStartAppending(tables[i].first.get()));
    }

    for (const auto& geometry : geometries) {
        geometry::BinaryBuffer wkb = geometry::to_wkb(geometry);
        geometry::Envelope envelope = geometry::envelope(geometry);

        ArrowBufferView wkb_view;
        wkb_view.data.data = wkb.data();
        wkb_view.size_bytes = static_cast<int64_t>(wkb.size());

        NANOARROW_THROW_NOT_OK(
            ArrowArrayAppendBytes(tables.front().first.get(), wkb_view));

        for (size_t i = 0; i < spatial_axes.size(); ++i) {
            NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(
                tables[2 * i + 1].first.get(), envelope.range.at(i).first));
            NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(
                tables[2 * i + 2].first.get(), envelope.range.at(i).second));
        }
    }

    for (size_t i = 0; i < tables.size(); ++i) {
        ArrowError error;
        NANOARROW_THROW_NOT_OK(
            ArrowArrayFinishBuildingDefault(tables[i].first.get(), &error));
    }

    return tables;
}

ArrowTable SOMAGeometryDataFrame::_reconstruct_geometry_data_table(
    ArrowTable original_data, const std::vector<ArrowTable>& wkb_data) {
    std::vector<std::string> spatial_axes = this->spatial_column_names();
    std::unordered_set<std::string> unique_column_names;
    std::unique_ptr<ArrowSchema> arrow_schema = std::make_unique<ArrowSchema>(
        ArrowSchema{});
    std::unique_ptr<ArrowArray> arrow_array = std::make_unique<ArrowArray>(
        ArrowArray{});

    for (int64_t i = 0; i < original_data.second->n_children; ++i) {
        unique_column_names.insert(original_data.second->children[i]->name);
    }
    for (size_t i = 0; i < wkb_data.size(); ++i) {
        unique_column_names.insert(wkb_data[i].second->name);
    }

    NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
        arrow_schema.get(), ArrowType::NANOARROW_TYPE_STRUCT));
    NANOARROW_THROW_NOT_OK(ArrowSchemaAllocateChildren(
        arrow_schema.get(), unique_column_names.size()));
    NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
        arrow_array.get(), ArrowType::NANOARROW_TYPE_STRUCT));
    NANOARROW_THROW_NOT_OK(ArrowArrayAllocateChildren(
        arrow_array.get(), unique_column_names.size()));

    // First add the wkb data columns so that already existing columns in the
    // original data except `soma_geometry` can overwrite the generated columns.

    for (size_t i = 0; i < wkb_data.size(); ++i) {
        ArrowSchemaMove(wkb_data[i].second.get(), arrow_schema->children[i]);
        ArrowArrayMove(wkb_data[i].first.get(), arrow_array->children[i]);
    }

    int64_t index = wkb_data.size();
    for (int64_t i = 0; i < original_data.second->n_children; ++i) {
        if (strcmp(original_data.second->children[i]->name, "soma_geometry") ==
            0) {
            continue;
        }

        bool replaced = false;
        for (size_t j = 0; j < wkb_data.size(); ++j) {
            if (strcmp(
                    arrow_schema->children[j]->name,
                    original_data.second->children[i]->name) == 0) {
                arrow_schema->children[j]->release(arrow_schema->children[j]);
                arrow_array->children[j]->release(arrow_array->children[j]);

                ArrowSchemaMove(
                    original_data.second->children[i],
                    arrow_schema->children[j]);
                ArrowArrayMove(
                    original_data.first->children[i], arrow_array->children[j]);

                replaced = true;
                break;
            }
        }

        if (!replaced) {
            ArrowSchemaMove(
                original_data.second->children[i],
                arrow_schema->children[index]);
            ArrowArrayMove(
                original_data.first->children[i], arrow_array->children[index]);

            ++index;
        }
    }

    return ArrowTable(std::move(arrow_array), std::move(arrow_schema));
}

}  // namespace tiledbsoma