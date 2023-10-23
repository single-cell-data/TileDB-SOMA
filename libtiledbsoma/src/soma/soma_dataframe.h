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

    SOMADataFrame(std::shared_ptr<SOMAArray> array)
        : array_(array){};

    SOMADataFrame() = delete;
    SOMADataFrame(const SOMADataFrame&) = default;
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
     * @brief Reset the state of this SOMADataFrame object to prepare for a
     * new query, while holding the array open.
     *
     * @param column_names
     * @param batch_size
     * @param result_order
     */
    void reset(
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        ResultOrder result_order = ResultOrder::automatic) {
        array_->reset(column_names, batch_size, result_order);
    }

    /**
     * @brief Check if the SOMADataFrame exists at the URI.
     */
    static bool exists(std::string_view uri);

    /**
     * Check if the SOMADataFrame is open.
     *
     * @return bool true if open
     */
    bool is_open() const;

    OpenMode mode() const {
        return array_->mode();
    }

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
     * Return optional timestamp pair SOMADataFrame was opened with.
     */
    std::optional<std::pair<uint64_t, uint64_t>> timestamp() {
        return array_->timestamp();
    }

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
    int64_t ndim() const;

    /**
     * Return the number of rows.
     *
     * @return int64_t
     */
    int64_t count() const;

    /**
     * @brief Get the capacity of each dimension.
     *
     * @return A vector with length equal to the number of dimensions; each
     * value in the vector is the capcity of each dimension.
     */
    std::vector<int64_t> shape() {
        return array_->shape();
    }

    /**
     * Retrieves the non-empty domain of the column index.
     *
     * @return int64_t
     */
    template <typename T>
    std::pair<T, T> non_empty_domain(const std::string& column_index_name) {
        return array_->non_empty_domain<T>(column_index_name);
    };

    /**
     * Retrieves the non-empty domain of the column index.
     * Applicable only to var-sized dimensions.
     */
    std::pair<std::string, std::string> non_empty_domain_var(
        const std::string& column_index_name) {
        return array_->non_empty_domain_var(column_index_name);
    };

    /**
     * Returns the domain of the given column index.
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    template <typename T>
    std::pair<T, T> domain(const std::string& column_index_name) const {
        return array_->domain<T>(column_index_name);
    }

    /**
     * @brief Read the next chunk of results from the query. If all results have
     * already been read, std::nullopt is returned.
     */
    std::optional<std::shared_ptr<ArrayBuffers>> read_next();

    /**
     * @brief Return true if `read_next` returned all results from the
     * query. The return value is false if the query was incomplete.
     *
     * @return True if last call to `read_next` returned all results of the
     * query
     */
    bool results_complete() {
        return array_->results_complete();
    }

    /**
     * @brief Set the dimension slice using one point
     *
     * @note Partitioning is not supported
     *
     * @tparam T
     * @param dim
     * @param point
     */
    template <typename T>
    void set_dim_point(const std::string& dim, const T& point) {
        array_->set_dim_point(dim, point);
    }

    /**
     * @brief Set the dimension slice using multiple points, with support
     * for partitioning.
     *
     * @tparam T
     * @param dim
     * @param points
     */
    template <typename T>
    void set_dim_points(
        const std::string& dim,
        const tcb::span<T> points,
        int partition_index,
        int partition_count) {
        array_->set_dim_points(dim, points, partition_index, partition_count);
    }

    /**
     * @brief Set the dimension slice using multiple points
     *
     * @note Partitioning is not supported
     *
     * @tparam T
     * @param dim
     * @param points
     */
    template <typename T>
    void set_dim_points(const std::string& dim, const std::vector<T>& points) {
        array_->set_dim_points(dim, points);
    }

    /**
     * @brief Set the dimension slice using multiple ranges
     *
     * @note Partitioning is not supported
     *
     * @tparam T
     * @param dim
     * @param ranges
     */
    template <typename T>
    void set_dim_ranges(
        const std::string& dim, const std::vector<std::pair<T, T>>& ranges) {
        array_->set_dim_ranges(dim, ranges);
    }

    /**
     * @brief Set a query condition.
     *
     * @param qc Query condition
     */
    void set_condition(QueryCondition& qc) {
        array_->set_condition(qc);
    }

    /**
     * @brief Returns the column names set by the query.
     *
     * @return std::vector<std::string>
     */
    std::vector<std::string> column_names() {
        return array_->column_names();
    }

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
