/**
 * @file   soma_point_cloud_dataframe.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAPointCloudDataFrame class.
 */

#ifndef SOMA_POINT_CLOUD_DATAFRAME
#define SOMA_POINT_CLOUD_DATAFRAME

#include <filesystem>

#include "soma_array.h"

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
     * @param ctx SOMAContext
     * @param platform_config Optional config parameter dictionary
     * @param timestamp Optional the timestamp range to write SOMA metadata info
     */
    static void create(
        std::string_view uri,
        std::unique_ptr<ArrowSchema> schema,
        ArrowTable index_columns,
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
        std::vector<std::string> column_names = {},
        ResultOrder result_order = ResultOrder::automatic,
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
     * @param column_names Columns to read
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp Timestamp
     */
    SOMAPointCloudDataFrame(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::vector<std::string> column_names,
        ResultOrder result_order,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMAArray(
              mode,
              uri,
              ctx,
              std::filesystem::path(uri).filename().string(),  // array name
              column_names,
              "auto",  // batch_size
              result_order,
              timestamp) {
    }

    SOMAPointCloudDataFrame(const SOMAArray& other)
        : SOMAArray(other) {
    }

    SOMAPointCloudDataFrame() = delete;
    SOMAPointCloudDataFrame(const SOMAPointCloudDataFrame&) = default;
    SOMAPointCloudDataFrame(SOMAPointCloudDataFrame&&) = delete;
    ~SOMAPointCloudDataFrame() = default;

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
};
}  // namespace tiledbsoma

#endif  // SOMA_POINT_CLOUD
