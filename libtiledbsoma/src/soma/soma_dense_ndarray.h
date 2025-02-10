/** * @file   soma_dense_ndarray.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMADenseNDArray class.
 */

#ifndef SOMA_DENSE_NDARRAY
#define SOMA_DENSE_NDARRAY

#include <filesystem>

#include "soma_array.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMADenseNDArray : public SOMAArray {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMADenseNDArray object at the given URI.
     *
     * @param uri URI to create the SOMASparseNDArray
     * @param format Arrow type to create the soma_data
     * @param index_columns The index column names with associated domains
     * and tile extents per dimension
     * @param ctx SOMAContext
     * @param platform_config Optional config parameter dictionary
     * @param timestamp Optional the timestamp range to write SOMA metadata info
     */
    static void create(
        std::string_view uri,
        std::string_view format,
        const ArrowTable& index_columns,
        std::shared_ptr<SOMAContext> ctx,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open and return a SOMADenseNDArray object at the given URI.
     *
     * @param uri URI to create the SOMADenseNDArray
     * @param mode read or write
     * @param ctx SOMAContext
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::unique_ptr<SOMADenseNDArray> SOMADenseNDArray
     */
    static std::unique_ptr<SOMADenseNDArray> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Check if the SOMADenseNDArray exists at the URI.
     *
     * @param uri URI to create the SOMADenseNDArray
     * @param ctx SOMAContext
     */
    static bool exists(std::string_view uri, std::shared_ptr<SOMAContext> ctx);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMADenseNDArray object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param timestamp Timestamp
     */
    SOMADenseNDArray(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp)
        : SOMAArray(mode, uri, ctx, timestamp) {
    }

    SOMADenseNDArray(const SOMAArray& other)
        : SOMAArray(other) {
    }

    SOMADenseNDArray() = delete;
    SOMADenseNDArray(const SOMADenseNDArray&) = default;
    SOMADenseNDArray(SOMADenseNDArray&&) = delete;
    ~SOMADenseNDArray() = default;

    using SOMAArray::open;

    /**
     * Return whether the SOMADenseNDArray is sparse.
     *
     * @return false
     */
    bool is_sparse() {
        return false;
    }

    /**
     * Return the data schema, in the form of an ArrowSchema.
     *
     * @return std::unique_ptr<ArrowSchema>
     */
    std::unique_ptr<ArrowSchema> schema() const;

    /**
     * @brief Get the soma_data's dtype in the form of an Arrow
     * format string.
     *
     * @return std::string_view Arrow format string.
     */
    std::string_view soma_data_type();
};
}  // namespace tiledbsoma

#endif  // SOMA_DENSE_NDARRAY
