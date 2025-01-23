/**
 * @file   soma_scene.h
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

#ifndef SOMA_SCENE
#define SOMA_SCENE

#include <tiledb/tiledb>

#include "soma_collection.h"
#include "soma_coordinates.h"

namespace tiledbsoma {

using namespace tiledb;
class SOMAScene : public SOMACollection {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAScene object at the given URI.
     *
     * @param uri URI to create the SOMAScene
     * @param schema TileDB ArraySchema
     * @param timestamp Optional pair indicating timestamp start and end
     */
    static void create(
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        const std::optional<SOMACoordinateSpace>& coordinate_space,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open a group at the specified URI and return SOMAScene
     * object.
     *
     * @param uri URI of the array
     * @param mode read or write
     * @param ctx TileDB context
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::shared_ptr<SOMAScene> SOMAScene
     */
    static std::unique_ptr<SOMAScene> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    SOMAScene(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMACollection(mode, uri, ctx, timestamp) {
    }

    SOMAScene(const SOMACollection& other)
        : SOMACollection(other) {
    }

    SOMAScene() = delete;
    SOMAScene(const SOMAScene&) = default;
    SOMAScene(SOMAScene&&) = default;
    ~SOMAScene() = default;

    inline const std::optional<SOMACoordinateSpace>& coordinate_space() const {
        return coord_space_;
    };

    /**
     * @brief Get the collection of imagery data.
     *
     * @return std::shared_ptr<SOMACollection>
     */
    std::shared_ptr<SOMACollection> img();

    /**
     * @brief Get the collection of observation location data.
     *
     * @return std::shared_ptr<SOMACollection>
     */
    std::shared_ptr<SOMACollection> obsl();

    /**
     * @brief Get the collection of collections of variable location data
     * for different measurements.
     *
     * @return std::shared_ptr<SOMACollection>
     */
    std::shared_ptr<SOMACollection> varl();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    std::optional<SOMACoordinateSpace> coord_space_ = std::nullopt;

    // A collection of imagery data.
    std::shared_ptr<SOMACollection> img_ = nullptr;

    // A collection of observation location data.
    std::shared_ptr<SOMACollection> obsl_ = nullptr;

    // A collection of collections of variable location data for measurements.
    std::shared_ptr<SOMACollection> varl_ = nullptr;
};
}  // namespace tiledbsoma

#endif  // SOMA_SCENE
