/**
 * @file   soma_sparse_ndarray.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMASparseNDArray class.
 */

#ifndef SOMA_SPARSE_NDARRAY
#define SOMA_SPARSE_NDARRAY

#include <filesystem>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include "../common/soma_column_selection.h"
#include "soma_array.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMASparseNDArray : public SOMAArray {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMASparseNDArray object at the given URI.
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
     * @brief Open and return a SOMASparseNDArray object at the given URI.
     *
     * @param uri URI to create the SOMASparseNDArray
     * @param mode read or write
     * @param ctx SOMAContext
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::unique_ptr<SOMASparseNDArray> SOMASparseNDArray
     */
    static std::unique_ptr<SOMASparseNDArray> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Check if the SOMASparseNDArray exists at the URI.
     *
     * @param ctx SOMAContext
     */
    static inline bool exists(std::string_view uri, std::shared_ptr<SOMAContext> ctx) {
        return SOMAArray::_exists(uri, "SOMASparseNDArray", ctx);
    }

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMASparseNDArray object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param timestamp Timestamp
     */
    SOMASparseNDArray(
        OpenMode mode, std::string_view uri, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp)
        : SOMAArray(mode, uri, ctx, timestamp) {
    }

    SOMASparseNDArray(const SOMAArray& other)
        : SOMAArray(other) {
    }

    SOMASparseNDArray() = delete;
    SOMASparseNDArray(const SOMASparseNDArray&) = default;
    SOMASparseNDArray(SOMASparseNDArray&&) = delete;
    ~SOMASparseNDArray() = default;

    using SOMAArray::open;

    /**
     * Return whether the SOMASparseNDArray is sparse.
     *
     * @return true
     */
    bool is_sparse() {
        return true;
    }

    /**
     * Return the data schema, in the form of an ArrowSchema.
     *
     * @return std::unique_ptr<ArrowSchema>
     */
    managed_unique_ptr<ArrowSchema> schema() const;

    /**
     * @brief Get the soma_data's dtype in the form of an Arrow
     * format string.
     *
     * @return std::string_view Arrow format string.
     */
    std::string_view soma_data_type();

    /**
     * @brief Delete cells from the array.
     *
     * This deletes cells that fall within the intersection of the provided coordinates. The coordinates are specified by
     * a per-dimension vector. If any coordinates are missing (length < number of dimensions) then the remaining dimensions
     * are set to be unconstrained. At least one constraint must be set.
     *
     * Acceptable ways to index a dimension:
     * * A slice that selects all coordinates from the starting value to the ending values including both end points. The slice
     *   must overlap the valid domain of the dimension.
     * * A list of points to select. All points must be contained inside the domain of the dimension.
     * * A default placeholder (`std::monostate`) - all values are selected (i.e. unconstrained).
     *
     * @param coords A per-dimension vector of coordinates - the intersection of which is the cells to delete.
     */
    void delete_cells(const std::vector<SOMAColumnSelection<int64_t>>& coords);
};
}  // namespace tiledbsoma

#endif  // SOMA_SPARSE_NDARRAY
