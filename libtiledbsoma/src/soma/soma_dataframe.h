/**
 * @file   soma_dataframe.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
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

#include <tiledb/tiledb>
#include "enums.h"
#include "soma_array.h"
#include "soma_object.h"

namespace tiledbsoma {

class SOMAArray;
class ArrayBuffers;

using namespace tiledb;

class SOMADataFrame : public SOMAObject {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMADataFrame object at the given URI.
     *
     * @param uri URI to create the SOMADataFrame
     * @param schema TileDB ArraySchema
     * @param platform_config Optional config parameter dictionary
     * @return std::shared_ptr<SOMADataFrame> opened in read mode
     */
    static std::unique_ptr<SOMADataFrame> create(
        std::string_view uri,
        ArraySchema schema,
        std::map<std::string, std::string> platform_config = {});

    /**
     * @brief Create a SOMADataFrame object at the given URI.
     *
     * @param uri URI to create the SOMADataFrame
     * @param schema TileDB ArraySchema
     * @param ctx TileDB context
     * @return std::shared_ptr<SOMADataFrame> opened in read mode
     */
    static std::unique_ptr<SOMADataFrame> create(
        std::string_view uri, ArraySchema schema, std::shared_ptr<Context> ctx);

    /**
     * @brief Open and return a SOMADataFrame object at the given URI.
     *
     * @param mode read or write
     * @param uri URI to create the SOMADataFrame
     * @param column_names A list of column names to use as user-defined index
     * columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must
     * exist in the schema, and at least one index column name is required.
     * @param platform_config Platform-specific options used to create this
     * DataFrame
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::shared_ptr<SOMADataFrame> SOMADataFrame
     */
    static std::unique_ptr<SOMADataFrame> open(
        std::string_view uri,
        OpenMode mode,
        std::map<std::string, std::string> platform_config = {},
        std::vector<std::string> column_names = {},
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Open and return a SOMADataFrame object at the given URI.
     *
     * @param mode read or write
     * @param ctx TileDB context
     * @param uri URI to create the SOMADataFrame
     * @param schema TileDB ArraySchema
     * @param column_names A list of column names to use as user-defined index
     * columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must
     * exist in the schema, and at least one index column name is required.
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::shared_ptr<SOMADataFrame> SOMADataFrame
     */
    static std::unique_ptr<SOMADataFrame> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names = {},
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

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
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        ResultOrder result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    SOMADataFrame() = delete;
    SOMADataFrame(const SOMADataFrame&) = delete;
    SOMADataFrame(SOMADataFrame&&) = default;
    ~SOMADataFrame() = default;

    /**
     * Open the SOMADataFrame object.
     *
     * @param mode read or write
     * @param timestamp Timestamp
     */
    void open(
        OpenMode mode,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * Close the SOMADataFrame object.
     */
    void close();

    /**
     * Check if the SOMADataFrame is open.
     *
     * @return bool true if open
     */
    bool is_open() const;

    /**
     * Return the constant "SOMADataFrame".
     *
     * @return std::string
     */
    const std::string type() const {
        return "SOMADataFrame";
    }

    /**
     * @brief Get the URI of the SOMADataFrame.
     *
     * @return std::string URI
     */
    const std::string uri() const;

    /**
     * Get the Context associated with the SOMADataFrame.
     *
     * @return std::shared_ptr<Context>
     */
    std::shared_ptr<Context> ctx();

    /**
     * Return the data schema, in the form of a TileDB ArraySchema.
     *
     * @return std::shared_ptr<ArraySchema>
     */
    std::shared_ptr<ArraySchema> schema() const;

    /**
     * Return the index (dimension) column names.
     *
     * @return std::vector<std::string>
     */
    const std::vector<std::string> index_column_names() const;

    /**
     * Return the number of rows in the SOMADataFrame.
     *
     * @return int64_t
     */
    int64_t count() const;

    /**
     * @brief Read the next chunk of results from the query. If all results have
     * already been read, std::nullopt is returned.
     */
    std::optional<std::shared_ptr<ArrayBuffers>> read_next();

    /**
     * @brief Write data to the dataframe.
     * @param buffers The ArrayBuffers to write
     */
    void write(std::shared_ptr<ArrayBuffers> buffers);

    /**
     * Set metadata key-value items to a SOMADataFrame. The SOMADataFrame must
     * opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be added. UTF-8 encodings
     *     are acceptable.
     * @param value_type The datatype of the value.
     * @param value_num The value may consist of more than one items of the
     *     same datatype. This argument indicates the number of items in the
     *     value component of the metadata.
     * @param value The metadata value in binary form.
     *
     * @note The writes will take effect only upon closing the array.
     */
    void set_metadata(
        const std::string& key,
        tiledb_datatype_t value_type,
        uint32_t value_num,
        const void* value) {
        array_->set_metadata(key, value_type, value_num, value);
    }

    /**
     * Delete a metadata key-value item from an open SOMADataFrame. The
     * SOMADataFrame must be opened in WRITE mode, otherwise the function will
     * error out.
     *
     * @param key The key of the metadata item to be deleted.
     *
     * @note The writes will take effect only upon closing the group.
     *
     * @note If the key does not exist, this will take no effect
     *     (i.e., the function will not error out).
     */
    void delete_metadata(const std::string& key) {
        array_->delete_metadata(key);
    }

    /**
     * @brief Given a key, get the associated value datatype, number of
     * values, and value in binary form.
     *
     * The value may consist of more than one items of the same datatype. Keys
     * that do not exist in the metadata will be return NULL for the value.
     *
     * **Example:**
     * @code{.cpp}
     * // Open the group for reading
     * tiledbsoma::SOMAGroup soma_group = SOMAGroup::open(TILEDB_READ,
     "s3://bucket-name/group-name");
     * tiledbsoma::MetadataValue meta_val = soma_group->get_metadata("key");
     * std::string key = std::get<MetadataInfo::key>(meta_val);
     * tiledb_datatype_t dtype = std::get<MetadataInfo::dtype>(meta_val);
     * uint32_t num = std::get<MetadataInfo::num>(meta_val);
     * const void* value = *((const
     int32_t*)std::get<MetadataInfo::value>(meta_val));
     * @endcode
     *
     * @param key The key of the metadata item to be retrieved. UTF-8 encodings
     *     are acceptable.
     * @return MetadataValue (std::tuple<std::string, tiledb_datatype_t,
     * uint32_t, const void*>)
     */
    std::optional<MetadataValue> get_metadata(const std::string& key) {
        return array_->get_metadata(key);
    }

    /**
     * Get a mapping of all metadata keys with its associated value datatype,
     * number of values, and value in binary form.
     *
     * @return std::map<std::string, MetadataValue>
     */
    std::map<std::string, MetadataValue> get_metadata() {
        return array_->get_metadata();
    }

    /**
     * Check if the key exists in metadata from an open SOMADataFrame.
     *
     * @param key The key of the metadata item to be checked. UTF-8 encodings
     *     are acceptable.
     * @return true if the key exists, else false.
     */
    bool has_metadata(const std::string& key) {
        return array_->has_metadata(key);
    }

    /**
     * Return then number of metadata items in an open SOMADataFrame. The group
     * must be opened in READ mode, otherwise the function will error out.
     */
    uint64_t metadata_num() const {
        return array_->metadata_num();
    }

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // SOMAArray
    std::shared_ptr<SOMAArray> array_;
};
}  // namespace tiledbsoma

#endif  // SOMA_DATAFRAME
