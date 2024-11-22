/**
 * @file   soma_geometry_dataframe.h
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
 *   This file defines the SOMAGeometryDataFrame class.
 */

#ifndef SOMA_GEOMETRY_DATAFRAME
#define SOMA_GEOMETRY_DATAFRAME

#include <filesystem>
#include <vector>

#include "soma_array.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMAGeometryDataFrame : virtual public SOMAArray {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAGeometryDataFrame object at the given URI.
     *
     * @param uri URI to create the SOMAGeometryDataFrame
     * @param schema Arrow schema
     * @param index_columns The index column names with associated domains
     * and tile extents per dimension
     * @param spatial_columns The spatial column names with associated domains
     * and tile extents per dimension
     * @param ctx SOMAContext
     * @param platform_config Optional config parameter dictionary
     * @param timestamp Optional the timestamp range to write SOMA metadata info
     */
    static void create(
        std::string_view uri,
        const std::unique_ptr<ArrowSchema>& schema,
        const ArrowTable& index_columns,
        const ArrowTable& spatial_columns,
        std::shared_ptr<SOMAContext> ctx,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open and return a SOMAGeometryDataFrame object at the given URI.
     *
     * @param uri URI to create the SOMAGeometryDataFrame
     * @param mode read or write
     * @param ctx SOMAContext
     * @param column_names A list of column names to use as user-defined index
     * columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must
     * exist in the schema, and at least one index column name is required.
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::unique_ptr<SOMAGeometryDataFrame> SOMAGeometryDataFrame
     */
    static std::unique_ptr<SOMAGeometryDataFrame> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::vector<std::string> column_names = {},
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Check if the SOMAGeometryDataFrame exists at the URI.
     *
     * @param uri URI to create the SOMAGeometryDataFrame
     * @param ctx SOMAContext
     */
    static bool exists(std::string_view uri, std::shared_ptr<SOMAContext> ctx);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMAGeometryDataFrame object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param column_names Columns to read
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp Timestamp
     */
    SOMAGeometryDataFrame(
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

    SOMAGeometryDataFrame(const SOMAArray& other)
        : SOMAArray(other) {
    }

    SOMAGeometryDataFrame() = delete;
    SOMAGeometryDataFrame(const SOMAGeometryDataFrame&) = default;
    SOMAGeometryDataFrame(SOMAGeometryDataFrame&&) = delete;
    ~SOMAGeometryDataFrame() = default;

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
     * Return the spatial column names.
     *
     * @return std::vector<std::string>
     */
    const std::vector<std::string> spatial_column_names() const;

    /**
     * Return the number of rows.
     *
     * @return int64_t
     */
    uint64_t count();

    void set_array_data(
        std::unique_ptr<ArrowSchema> arrow_schema,
        std::unique_ptr<ArrowArray> arrow_array) override;

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    /**
     * @brief Cast an array containing the outer rings of polygons to an Arrow
     * array holding the WKB encoded polygons and generate the additional index
     * column arrays based on the spatial axes.
     */
    std::vector<ArrowTable> _cast_polygon_vertex_list_to_wkb(ArrowArray* array);

    /**
     * @brief Create a new ArrowTable by merging the generated WKB and spatial
     * index arrays and the original data.
     *
     * @remark Generated columns have predefined names. Any generated column
     * with name already present in the original data will be skipped.
     */
    ArrowTable _reconstruct_geometry_data_table(
        ArrowTable original_data, const std::vector<ArrowTable>& wkb_data);
};
}  // namespace tiledbsoma

#endif  // SOMA_GEOMETRY_DATAFRAME
