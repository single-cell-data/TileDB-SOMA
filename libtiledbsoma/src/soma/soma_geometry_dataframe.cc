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
#include "../utils/transformer.h"
#include "../utils/util.h"
#include "soma_geometry_column.h"
#include "soma_transformers.h"

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
    std::optional<TimestampRange> timestamp) {
    return std::make_unique<SOMAGeometryDataFrame>(mode, uri, ctx, timestamp);
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

}  // namespace tiledbsoma
