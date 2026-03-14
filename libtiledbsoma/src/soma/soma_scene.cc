/**
 * @file   soma_scene.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAScene class.
 */

#include "soma_scene.h"
#include "soma_collection.h"
#include "utils/common.h"

#include "common/datatype/datatype.h"
#include "common/datatype/utils.h"

namespace tiledbsoma {
using namespace common::type;

//===================================================================
//= public static
//===================================================================

void SOMAScene::create(
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    const std::optional<SOMACoordinateSpace>& coordinate_space,
    std::optional<TimestampRange> timestamp) {
    std::unordered_map<std::string, std::string> schema_metadata{
        {SPATIAL_ENCODING_VERSION_KEY, SPATIAL_ENCODING_VERSION_VAL}};
    if (coordinate_space.has_value()) {
        schema_metadata[SOMA_COORDINATE_SPACE_KEY] = coordinate_space->to_string();
    }
    try {
        SOMAGroup::create(ctx, uri, "SOMAScene", schema_metadata, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAScene> SOMAScene::open(
    std::string_view uri, OpenMode mode, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    try {
        return std::make_unique<SOMAScene>(mode, uri, ctx, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

SOMAScene::SOMAScene(
    OpenMode mode, std::string_view uri, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp)
    : SOMACollectionBase(mode, uri, ctx, timestamp, "SOMAScene") {
    auto coord_space_metadata = get_metadata(SOMA_COORDINATE_SPACE_KEY);
    if (coord_space_metadata.has_value()) {
        coord_space_ = SOMACoordinateSpace::from_metadata(coord_space_metadata.value());
    }
}

std::shared_ptr<SOMACollection> SOMAScene::img() {
    return std::dynamic_pointer_cast<SOMACollection>(get("img"));
}

std::shared_ptr<SOMACollection> SOMAScene::obsl() {
    return std::dynamic_pointer_cast<SOMACollection>(get("obsl"));
}

std::shared_ptr<SOMACollection> SOMAScene::varl() {
    return std::dynamic_pointer_cast<SOMACollection>(get("varl"));
}

std::string SOMAScene::classname() const {
    return "Scene";
}

}  // namespace tiledbsoma
