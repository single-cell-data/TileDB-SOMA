/**
 * @file   soma_geometry_dataframe.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
#include "soma_geometry_column.h"

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
    const SOMACoordinateSpace& coordinate_space,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    auto [tiledb_schema, soma_schema_extension] =
        ArrowAdapter::tiledb_schema_from_arrow_schema(
            ctx->tiledb_ctx(),
            schema,
            index_columns,
            std::make_optional(coordinate_space),
            "SOMAGeometryDataFrame",
            true,
            platform_config);

    auto array = SOMAArray::_create(
        ctx,
        uri,
        tiledb_schema,
        "SOMAGeometryDataFrame",
        soma_schema_extension.dump(),
        timestamp);

    // Add additional geometry dataframe metadata.
    array.put_metadata(
        SPATIAL_ENCODING_VERSION_KEY,
        TILEDB_STRING_UTF8,
        static_cast<uint32_t>(SPATIAL_ENCODING_VERSION_VAL.size()),
        SPATIAL_ENCODING_VERSION_VAL.c_str());
    const auto coord_space_metadata = coordinate_space.to_string();
    array.put_metadata(
        SOMA_COORDINATE_SPACE_KEY,
        TILEDB_STRING_UTF8,
        static_cast<uint32_t>(coord_space_metadata.size()),
        coord_space_metadata.c_str());
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

uint64_t SOMAGeometryDataFrame::count() {
    return this->nnz();
}

void SOMAGeometryDataFrame::set_array_data(
    std::unique_ptr<ArrowSchema> arrow_schema,
    std::unique_ptr<ArrowArray> arrow_array) {
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

void SOMAGeometryDataFrame::initialize() {
    auto coordinate_space_meta = get_metadata(SOMA_COORDINATE_SPACE_KEY);

    if (!coordinate_space_meta.has_value()) {
        throw TileDBSOMAError(std::format(
            "[SOMAGeometryDataFrame][initialize] Missing required '{}' "
            "metadata key.",
            SOMA_COORDINATE_SPACE_KEY));
    }

    coord_space_ = std::apply(
        SOMACoordinateSpace::from_metadata, coordinate_space_meta.value());
}

std::vector<ArrowTable> SOMAGeometryDataFrame::_cast_polygon_vertex_list_to_wkb(
    ArrowArray* array) {
    // Initialize a vector to hold all the Arrow tables containing the
    // transformed geometry data
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

    for (size_t i = 0; i < coord_space_.size(); ++i) {
        const auto axis = coord_space_.axis(i);

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
            (SOMA_GEOMETRY_DIMENSION_PREFIX + axis.name + "__min").c_str()));

        // Max spatial axis
        tables.push_back(ArrowTable(
            std::make_unique<ArrowArray>(), std::make_unique<ArrowSchema>()));
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            tables.back().first.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
            tables.back().second.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaSetName(
            tables.back().second.get(),
            (SOMA_GEOMETRY_DIMENSION_PREFIX + axis.name + "__max").c_str()));
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

        for (size_t i = 0; i < coord_space_.size(); ++i) {
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
