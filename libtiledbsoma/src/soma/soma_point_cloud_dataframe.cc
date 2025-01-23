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

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAPointCloudDataFrame::create(
    std::string_view uri,
    const std::unique_ptr<ArrowSchema>& schema,
    const ArrowTable& index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    auto tiledb_schema = ArrowAdapter::tiledb_schema_from_arrow_schema(
        ctx->tiledb_ctx(),
        schema,
        index_columns,
        "SOMAPointCloudDataFrame",
        true,
        platform_config);
    SOMAArray::create(
        ctx, uri, tiledb_schema, "SOMAPointCloudDataFrame", timestamp);
}

std::unique_ptr<SOMAPointCloudDataFrame> SOMAPointCloudDataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp) {
    auto array = std::make_unique<SOMAPointCloudDataFrame>(
        mode, uri, ctx, column_names, result_order, timestamp);

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
