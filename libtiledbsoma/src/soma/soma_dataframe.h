/**
 * @file   soma_dataframe.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023-2024 TileDB, Inc.
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
        std::unique_ptr<ArrowSchema> schema,
        ArrowTable index_columns,
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
     * @param column_names A list of column names to use as user-defined index
     * columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must
     * exist in the schema, and at least one index column name is required.
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::unique_ptr<SOMADataFrame> SOMADataFrame
     */
    static std::unique_ptr<SOMADataFrame> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::vector<std::string> column_names = {},
        ResultOrder result_order = ResultOrder::automatic,
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
        std::string_view name,
        std::map<std::string, std::string> platform_config);

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
     * o add_enmrs: attr_name -> pair of:
     *   - Arrow type string the value type, e.g. "f" or "U"
     *   - bool ordered
     */
    static void update_dataframe_schema(
        std::string uri,
        std::shared_ptr<SOMAContext> ctx,
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
     * @param column_names Columns to read
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp Timestamp
     */
    SOMADataFrame(
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
     * @param name Name of the array
     * @param platform_config Config parameter dictionary
     * @return std::unique_ptr<SOMADataFrame> SOMADataFrame
     */
    SOMADataFrame(
        OpenMode mode,
        std::string_view uri,
        std::string_view name,
        std::map<std::string, std::string> platform_config)
        : SOMAArray(
              mode,
              uri,
              std::make_shared<SOMAContext>(platform_config),
              name,
              {},
              "auto",
              ResultOrder::automatic) {
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
