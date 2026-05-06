/**
 * @file   soma_collection_base.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMACollectionBase class.
 */

#ifndef SOMA_COLLECTION_BASE
#define SOMA_COLLECTION_BASE

#include <tiledb/tiledb>

#include "../tiledb_adapter/platform_config.h"
#include "common/arrow/utils.h"
#include "enums.h"
#include "soma_group.h"
#include "utils/common.h"

namespace tiledbsoma {

class SOMACollection;
class SOMADataFrame;
class SOMADenseNDArray;
class SOMAExperiment;
class SOMAMeasurement;
class SOMASparseNDArray;

using namespace tiledb;

class SOMACollectionBase : public SOMAGroup {
   public:
    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMACollectionBase object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param key key of the array
     * @param timestamp Optional pair indicating timestamp start and end
     */
    SOMACollectionBase(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp,
        std::optional<std::string> soma_type);

    SOMACollectionBase(SOMAGroup&& other);

    SOMACollectionBase() = delete;
    SOMACollectionBase(const SOMACollectionBase&) = default;
    SOMACollectionBase(SOMACollectionBase&&) = default;
    virtual ~SOMACollectionBase();

    /**
     * Open the SOMAGroup object.
     *
     * @param mode read or write
     * @param timestamp Optional pair indicating timestamp start and end
     */
    void open(OpenMode mode, std::optional<TimestampRange> timestamp = std::nullopt) override;

    /**
     * Open the SOMAGroup object.
     *
     * @param mode read or write
     * @param timestamp Optional pair indicating timestamp start and end
     */
    void reopen(OpenMode mode, std::optional<TimestampRange> timestamp = std::nullopt) override;

    /**
     * Closes the SOMACollectionBase object.
     */
    void close([[maybe_unused]] bool recursive = false) override;

    /**
     * Get the SOMAObject associated with the key.
     *
     * @param key of member
     */
    std::shared_ptr<SOMAObject> get(const std::string& key);

    void set(
        const std::string& uri,
        URIType uri_type,
        const std::string& name,
        const std::string& soma_type,
        std::shared_ptr<SOMAObject> member,
        bool managed = true);

    /**
     * Remove a named member from the SOMAGroup.
     *
     * @param name of member
     */
    void del(const std::string& name) override;

    /**
     * Create and add a SOMACollection to the SOMACollectionBase.
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
     * Create and add a SOMAExperiment to the SOMACollectionBase.
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
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMAMeasurement to the SOMACollectionBase.
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
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMADataFrame to the SOMACollectionBase.
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
        const common::arrow::managed_unique_ptr<ArrowSchema>& schema,
        const common::arrow::ArrowTable& index_columns,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMADenseNDArray to the SOMACollectionBase.
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
        const common::arrow::ArrowTable& index_columns,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Create and add a SOMASparseNDArray to the SOMACollectionBase.
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
        const common::arrow::ArrowTable& index_columns,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    std::ostream& print(
        std::ostream& stream, int level = 0, std::optional<std::string> key = std::nullopt) const override;

   protected:
    //===================================================================
    //= protected non-static
    //===================================================================

    // Members of the SOMACollectionBase
    std::map<std::string, std::shared_ptr<SOMAObject>> children_;

   private:
    std::map<std::string, std::shared_ptr<std::once_flag>> flags_;
    std::unordered_set<std::string> managed_children_;
};
}  // namespace tiledbsoma

#endif  // SOMA_COLLECTION
