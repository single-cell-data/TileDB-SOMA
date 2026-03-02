/**
 * @file   soma_multiscale_image.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAMultiscaleImage class.
 */

#include "soma_multiscale_image.h"
#include "soma_collection.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAMultiscaleImage::create(
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    const SOMACoordinateSpace& coordinate_space,
    std::optional<TimestampRange> timestamp) {
    try {
        SOMAGroup::create(
            ctx,
            uri,
            "SOMAMultiscaleImage",
            {{SPATIAL_ENCODING_VERSION_KEY, SPATIAL_ENCODING_VERSION_VAL},
             {SOMA_COORDINATE_SPACE_KEY, coordinate_space.to_string()}},
            timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAMultiscaleImage> SOMAMultiscaleImage::open(
    std::string_view uri, OpenMode mode, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    try {
        return std::make_unique<SOMAMultiscaleImage>(mode, uri, ctx, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

}  // namespace tiledbsoma
