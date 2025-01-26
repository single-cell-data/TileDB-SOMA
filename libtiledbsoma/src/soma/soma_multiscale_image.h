/**
 * @file   soma_multiscale_image.h
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

#ifndef SOMA_MULTISCALE_IMAGE
#define SOMA_MULTISCALE_IMAGE

#include <tiledb/tiledb>

#include "soma_collection.h"
#include "soma_coordinates.h"

namespace tiledbsoma {

using namespace tiledb;
class SOMAMultiscaleImage : public SOMACollection {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAMultiscaleImage object at the given URI.
     *
     * @param uri URI to create the SOMAMultiscaleImage
     * @param schema TileDB ArraySchema
     * @param platform_config Optional config parameter dictionary
     */
    static void create(
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        const SOMACoordinateSpace& coordinate_space,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open a group at the specified URI and return SOMAMultiscaleImage
     * object.
     *
     * @param uri URI of the array
     * @param mode read or write
     * @param ctx TileDB context
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::shared_ptr<SOMAMultiscaleImage> SOMAMultiscaleImage
     */
    static std::unique_ptr<SOMAMultiscaleImage> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    SOMAMultiscaleImage(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMACollection(mode, uri, ctx, timestamp) {
    }

    SOMAMultiscaleImage(const SOMACollection& other)
        : SOMACollection(other) {
    }

    SOMAMultiscaleImage() = delete;
    SOMAMultiscaleImage(const SOMAMultiscaleImage&) = default;
    SOMAMultiscaleImage(SOMAMultiscaleImage&&) = default;
    ~SOMAMultiscaleImage() = default;

    inline const SOMACoordinateSpace& coordinate_space() const {
        return coord_space_;
    }

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    SOMACoordinateSpace coord_space_;
};
}  // namespace tiledbsoma

#endif  // SOMA_MULTISCALE_IMAGE
