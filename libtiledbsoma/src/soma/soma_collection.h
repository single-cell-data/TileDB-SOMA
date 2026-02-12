/**
 * @file   soma_collection.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMACollection class.
 */

#ifndef SOMA_COLLECTION
#define SOMA_COLLECTION

#include "../tiledb_adapter/platform_config.h"
#include "enums.h"
#include "soma_collection_base.h"
#include "utils/common.h"

namespace tiledbsoma {

class SOMACollection : public SOMACollectionBase {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMACollection object at the given URI.
     *
     * @param ctx TileDB context
     * @param uri URI to create the SOMACollection
     */
    static void create(
        std::string_view uri, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open a group at the specified URI and return SOMACollection
     * object.
     *
     * @param uri URI of the array
     * @param mode read or write
     * @param ctx TileDB context
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::shared_ptr<SOMACollection> SOMACollection
     */
    static std::unique_ptr<SOMACollection> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    using SOMACollectionBase::open;

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMACollection object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param key key of the array
     * @param timestamp Optional pair indicating timestamp start and end
     */
    SOMACollection(
        OpenMode mode, std::string_view uri, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp)
        : SOMACollectionBase(mode, uri, ctx, timestamp, "SOMACollection") {
    }

    SOMACollection(const SOMACollectionBase& other)
        : SOMACollectionBase(other) {
    }

    SOMACollection() = delete;
    SOMACollection(const SOMACollection&) = default;
    SOMACollection(SOMACollection&&) = default;
    virtual ~SOMACollection() = default;
};
}  // namespace tiledbsoma

#endif  // SOMA_COLLECTION
