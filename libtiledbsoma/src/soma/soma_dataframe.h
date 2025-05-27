/**
 * @file   soma_dataframe.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMADataFrame class.
 */

#ifndef SOMA_DATAFRAME
#define SOMA_DATAFRAME

#include <filesystem>

#include "soma_array.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMADataFrame : public SOMAArray {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMADataFrame object at the given URI.
     *
     * @param uri URI to create the SOMADataFrame
     * @param schema Arrow schema
     * @param index_columns The index column names with associated domains
     * and tile extents per dimension
     * @param ctx SOMAContext
     * @param platform_config Optional config parameter dictionary
     * @param timestamp Optional the timestamp range to write SOMA metadata info
     */
    static void create(
        std::string_view uri,
        const std::unique_ptr<ArrowSchema>& schema,
        const ArrowTable& index_columns,
        std::shared_ptr<SOMAContext> ctx,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open and return a SOMADataFrame object at the given URI.
     *
     * Note: the first two arguments uri and mode are reversed from
     * the SOMAArrayConstructor. This is an intentional decision to
     * avoid ambiguous-overload compiler errors. Even though
     * SOMADataFrame extends SOMAArray, callers using open
     * and wishing to obtain a SOMADataFrame rather than a SOMAArray
     * are advised to place the uri argument before the mode argument.
     *
     * @param uri URI to create the SOMADataFrame
     * @param mode read or write
     * @param ctx SOMAContext
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::unique_ptr<SOMADataFrame> SOMADataFrame
     */
    static std::unique_ptr<SOMADataFrame> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open and return a SOMADataFrame object at the given URI.
     *
     * This is nominally for TileDB-SOMA R use, since while we have
     * TileDB-SOMA-R and TileDB-R co-existing, it's less desirable
     * to pass ctx from one copy of core to another, and more
     * desirable to pass a config map.
     *
     * @param uri URI to create the SOMADataFrame
     * @param mode read or write
     * @param name Name of the array
     * @param platform_config Config parameter dictionary
     * @return std::unique_ptr<SOMADataFrame> SOMADataFrame
     */
    static std::unique_ptr<SOMADataFrame> open(
        std::string_view uri,
        OpenMode mode,
        std::map<std::string, std::string> platform_config,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Check if the SOMADataFrame exists at the URI.
     *
     * @param uri URI to create the SOMADataFrame
     * @param ctx SOMAContext
     */
    static bool exists(std::string_view uri, std::shared_ptr<SOMAContext> ctx);

    /**
     * This is for schema evolution.
     *
     * For non-enum attrs:
     *
     * o drop_cols: attr_name
     * o add_attrs: attr_name -> Arrow type string like "i" or "U"
     * o add_enmrs: no key present
     *
     * Enum attrs:
     *
     * o drop_cols: attr_name
     * o add_attrs: attr_name -> Arrow type string for the index
     *   type, e.g. 'c' for int8
     * o add_enmrs: attr_name -> tuple of:
     *   - Arrow type string the value type, e.g. "f" or "U"
     *   - bool ordered
     */
    void update_dataframe_schema(
        std::vector<std::string> drop_attrs,
        std::map<std::string, std::string> add_attrs,
        std::map<std::string, std::pair<std::string, bool>> add_enmrs);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMADataFrame object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param timestamp Timestamp
     */
    SOMADataFrame(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMAArray(mode, uri, ctx, timestamp) {
    }

    /**
     * @brief Construct a new SOMADataFrame object.
     *
     * This is nominally for TileDB-SOMA R use, since while we have
     * TileDB-SOMA-R and TileDB-R co-existing, it's less desirable
     * to pass ctx from one copy of core to another, and more
     * desirable to pass a config map.
     *
     * @param uri URI to create the SOMADataFrame
     * @param mode read or write
     * @param platform_config Config parameter dictionary
     * @return std::unique_ptr<SOMADataFrame> SOMADataFrame
     */
    SOMADataFrame(
        OpenMode mode,
        std::string_view uri,
        std::map<std::string, std::string> platform_config,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMAArray(
              mode,
              uri,
              std::make_shared<SOMAContext>(platform_config),
              timestamp) {
    }

    SOMADataFrame(const SOMAArray& other)
        : SOMAArray(other) {
    }

    SOMADataFrame() = delete;
    SOMADataFrame(const SOMADataFrame&) = default;
    SOMADataFrame(SOMADataFrame&&) = delete;
    ~SOMADataFrame() = default;

    using SOMAArray::open;

    /**
     * Return the data schema, in the form of a ArrowSchema.
     *
     * @return std::unique_ptr<ArrowSchema>
     */
    std::unique_ptr<ArrowSchema> schema() const;

    /**
     * Return the index (dimension) column names.
     *
     * @return std::vector<std::string>
     */
    const std::vector<std::string> index_column_names() const;

    /**
     * Return the number of rows.
     *
     * @return int64_t
     */
    uint64_t count();

    /**
     * While application-level SOMA DataFrame doesn't have shape
     * and maxshape, these are important test-point accessors,
     * as well as crucial for experiment-level resize within tiledbsoma.io.
     *
     * Note that the SOMA spec for SOMADataFrame mandates a .domain() accessor,
     * which is distinct, and type-polymorphic.
     *
     * @return std::optional<int64_t>
     */
    std::optional<int64_t> maybe_soma_joinid_shape();

    /**
     * See comments for maybe_soma_joinid_shape.
     *
     * @return std::optional<int64_t>
     */
    std::optional<int64_t> maybe_soma_joinid_maxshape();
};

}  // namespace tiledbsoma

#endif  // SOMA_DATAFRAME
