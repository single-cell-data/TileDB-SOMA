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

#include <tiledb/tiledb>

#include "enums.h"
#include "soma_dataframe.h"
#include "soma_dense_ndarray.h"
#include "soma_group.h"
#include "soma_object.h"
#include "soma_sparse_ndarray.h"

namespace tiledbsoma {

class SOMAExperiment;
class SOMAMeasurement;

using namespace tiledb;

class SOMACollection : public SOMAGroup {
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
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

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
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp)
        : SOMAGroup(
              mode,
              uri,
              ctx,
              std::filesystem::path(uri).filename().string(),  // group name
              timestamp){};

    SOMACollection(const SOMAGroup& other)
        : SOMAGroup(other) {
    }

    SOMACollection() = delete;
    SOMACollection(const SOMACollection&) = default;
    SOMACollection(SOMACollection&&) = default;
    virtual ~SOMACollection() = default;

    using iterator =
        typename std::map<std::string, std::shared_ptr<SOMAObject>>::iterator;
    iterator begin() {
        return children_.begin();
    }
    iterator end() {
        return children_.end();
    }

    using SOMAGroup::open;

    /**
     * Closes the SOMACollection object.
     */
    void close();

    /**
     * Get the SOMAObject associated with the key.
     *
     * @param key of member
     */
    std::unique_ptr<SOMAObject> get(const std::string& key);

    /**
     * Create and add a SOMACollection to the SOMACollection.
     *
     * @param key of collection
     * @param uri of SOMACollection to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param timestamp Optional the timestamp range to open at
     */
    std::shared_ptr<SOMACollection> add_new_collection(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMAExperiment to the SOMACollection.
     *
     * @param key of collection
     * @param uri of SOMAExperiment to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param timestamp Optional the timestamp range to open at
     */
    std::shared_ptr<SOMAExperiment> add_new_experiment(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        std::unique_ptr<ArrowSchema> schema,
        ArrowTable index_columns,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMAMeasurement to the SOMACollection.
     *
     * @param key of collection
     * @param uri of SOMAMeasurement to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param timestamp Optional the timestamp range to open at
     */
    std::shared_ptr<SOMAMeasurement> add_new_measurement(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        std::unique_ptr<ArrowSchema> schema,
        ArrowTable index_columns,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMADataFrame to the SOMACollection.
     *
     * @param key of dataframe
     * @param uri of SOMADataFrame to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param ctx SOMAContext
     * @param format Arrow type to create the soma_data
     * @param index_columns The index column names with associated domains
     * and tile extents per dimension
     * @param platform_config Optional config parameter dictionary
     * @param timestamp Optional the timestamp range to open at
     */
    std::shared_ptr<SOMADataFrame> add_new_dataframe(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        std::unique_ptr<ArrowSchema> schema,
        ArrowTable index_columns,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMADenseNDArray to the SOMACollection.
     *
     * @param key of dense array
     * @param uri of SOMADenseNDArray to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param ctx SOMAContext
     * @param format Arrow type to create the soma_data
     * @param index_columns The index column names with associated domains
     * and tile extents per dimension
     * @param platform_config Optional config parameter dictionary
     * @param timestamp Optional the timestamp range to open at
     */
    std::shared_ptr<SOMADenseNDArray> add_new_dense_ndarray(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        std::string_view format,
        ArrowTable index_columns,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMASparseNDArray to the SOMACollection.
     *
     * @param key of sparse array
     * @param uri of SOMASparseNDArray to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param ctx SOMAContext
     * @param format Arrow type to create the soma_data
     * @param index_columns The index column names with associated domains
     * and tile extents per dimension
     * @param platform_config Optional config parameter dictionary
     * @param timestamp Optional the timestamp range to open at
     */
    std::shared_ptr<SOMASparseNDArray> add_new_sparse_ndarray(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        std::string_view format,
        ArrowTable index_columns,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

   protected:
    //===================================================================
    //= protected non-static
    //===================================================================

    // Members of the SOMACollection
    std::map<std::string, std::shared_ptr<SOMAObject>> children_;
};
}  // namespace tiledbsoma

#endif  // SOMA_COLLECTION
