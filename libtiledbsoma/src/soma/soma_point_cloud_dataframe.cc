/**
 * @file   soma_point_cloud_dataframe.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAPointCloudDataFrame class.
 */

#include "soma_point_cloud_dataframe.h"
#include <tiledb/tiledb>
#include "utils/common.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAPointCloudDataFrame::create(
    std::string_view uri,
    const std::unique_ptr<ArrowSchema>& schema,
    const ArrowTable& index_columns,
    const SOMACoordinateSpace& coordinate_space,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    // Create TileDB array that is open for writing.
    auto [tiledb_schema, soma_schema_extension] =
        ArrowAdapter::tiledb_schema_from_arrow_schema(
            ctx->tiledb_ctx(),
            schema,
            index_columns,
            std::make_optional(coordinate_space),
            "SOMAPointCloudDataFrame",
            true,
            platform_config);
    auto array = SOMAArray::_create(
        ctx,
        uri,
        tiledb_schema,
        "SOMAPointCloudDataFrame",
        std::nullopt,
        timestamp);

    // Add additional point cloud dataframe metadata.
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

std::unique_ptr<SOMAPointCloudDataFrame> SOMAPointCloudDataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    auto array = std::make_unique<SOMAPointCloudDataFrame>(
        mode, uri, ctx, timestamp);

    if (!array->check_type("SOMAPointCloudDataFrame")) {
        throw TileDBSOMAError(
            "[SOMAPointCloudDataFrame::open] Object is not a "
            "SOMAPointCloudDataFrame");
    }

    return array;
}

bool SOMAPointCloudDataFrame::exists(
    std::string_view uri, std::shared_ptr<SOMAContext> ctx) {
    try {
        auto obj = SOMAObject::open(uri, OpenMode::read, ctx);
        return "SOMAPointCloudDataFrame" == obj->type();
    } catch (TileDBSOMAError& e) {
        return false;
    }
}

//===================================================================
//= public non-static
//===================================================================

std::unique_ptr<ArrowSchema> SOMAPointCloudDataFrame::schema() const {
    return this->arrow_schema();
}

const std::vector<std::string> SOMAPointCloudDataFrame::index_column_names()
    const {
    return this->dimension_names();
}

uint64_t SOMAPointCloudDataFrame::count() {
    return this->nnz();
}

}  // namespace tiledbsoma
