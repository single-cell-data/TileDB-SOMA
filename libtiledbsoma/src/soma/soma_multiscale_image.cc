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
        std::filesystem::path image_uri(uri);
        auto group = SOMAGroup::create(
            ctx, image_uri.string(), "SOMAMultiscaleImage", timestamp);

        // Set spatial-encoding metadata.
        group->set_metadata(
            SPATIAL_ENCODING_VERSION_KEY,
            TILEDB_STRING_UTF8,
            static_cast<uint32_t>(SPATIAL_ENCODING_VERSION_VAL.size()),
            SPATIAL_ENCODING_VERSION_VAL.c_str(),
            true);

        // Set coordinate space metadata.
        const auto coord_space_metadata = coordinate_space.to_string();
        group->set_metadata(
            SOMA_COORDINATE_SPACE_KEY,
            TILEDB_STRING_UTF8,
            static_cast<uint32_t>(coord_space_metadata.size()),
            coord_space_metadata.c_str(),
            true);

    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAMultiscaleImage> SOMAMultiscaleImage::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    try {
        auto group = std::make_unique<SOMAMultiscaleImage>(
            mode, uri, ctx, timestamp);

        if (!group->check_type("SOMAMultiscaleImage")) {
            throw TileDBSOMAError(
                "[SOMAMultiscaleImage::open] Object is not a "
                "SOMAMultiscaleImage");
        }

        return group;
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

}  // namespace tiledbsoma
