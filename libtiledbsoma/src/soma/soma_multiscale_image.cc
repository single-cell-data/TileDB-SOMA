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
    std::optional<TimestampRange> timestamp) {
    try {
        std::filesystem::path image_uri(uri);
        SOMAGroup::create(
            ctx, image_uri.string(), "SOMAMultiscaleImage", timestamp);
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
