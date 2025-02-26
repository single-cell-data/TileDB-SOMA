/**
 * @file   soma_point_cloud_dataframe.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAPointCloudDataFrame class.
 */

#ifndef SOMA_POINT_CLOUD_DATAFRAME
#define SOMA_POINT_CLOUD_DATAFRAME

#include <filesystem>

#include "soma_array.h"
#include "soma_coordinates.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMAPointCloudDataFrame : public SOMAArray {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAPointCloudDataFrame object at the given URI.
     *
     * @param uri URI to create the SOMAPointCloudDataFrame
     * @param schema Arrow schema
     * @param index_columns The index column names with associated domains
     * and tile extents per dimension
     * @param coordinate_space The coordinate space the PointCloudDataFrame
     * spatial axes are defined on.
     * @param ctx SOMAContext
     * @param platform_config Optional config parameter dictionary
     * @param timestamp Optional the timestamp range to write SOMA metadata info
     */
    static void create(
        std::string_view uri,
        const std::unique_ptr<ArrowSchema>& schema,
        const ArrowTable& index_columns,
        const SOMACoordinateSpace& coordinate_space,
        std::shared_ptr<SOMAContext> ctx,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open and return a SOMAPointCloudDataFrame object at the given URI.
     *
     * @param uri URI to create the SOMAPointCloudDataFrame
     * @param mode read or write
     * @param ctx SOMAContext
     * @param column_names A list of column names to use as user-defined index
     * columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must
     * exist in the schema, and at least one index column name is required.
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::unique_ptr<SOMAPointCloudDataFrame> SOMAPointCloudDataFrame
     */
    static std::unique_ptr<SOMAPointCloudDataFrame> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Check if the SOMAPointCloudDataFrame exists at the URI.
     *
     * @param uri URI to create the SOMAPointCloudDataFrame
     * @param ctx SOMAContext
     */
    static bool exists(std::string_view uri, std::shared_ptr<SOMAContext> ctx);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMAPointCloudDataFrame object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param timestamp Timestamp
     */
    SOMAPointCloudDataFrame(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMAArray(mode, uri, ctx, timestamp) {
    }

    SOMAPointCloudDataFrame(const SOMAArray& other)
        : SOMAArray(other) {
    }

    SOMAPointCloudDataFrame() = delete;
    SOMAPointCloudDataFrame(const SOMAPointCloudDataFrame&) = default;
    SOMAPointCloudDataFrame(SOMAPointCloudDataFrame&&) = delete;
    ~SOMAPointCloudDataFrame() = default;

    using SOMAArray::open;

    inline const SOMACoordinateSpace& coordinate_space() const {
        return coord_space_;
    };

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

   private:
    SOMACoordinateSpace coord_space_;
};
}  // namespace tiledbsoma

#endif  // SOMA_POINT_CLOUD
