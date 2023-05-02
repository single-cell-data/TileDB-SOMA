/**
 * @file   soma_array.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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
 *   This declares the SOMAArray
 */

#ifndef SOMA_ARRAY
#define SOMA_ARRAY

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <future>

#include <tiledb/tiledb>

#include "managed_query.h"

namespace tiledbsoma {
using namespace tiledb;

using MetadataValue =
    std::tuple<std::string, tiledb_datatype_t, uint32_t, const void*>;
enum MetadataInfo { key = 0, dtype, num, value };

class SOMAArray {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open an array at the specified URI and return SOMAArray
     * object.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param uri URI of the array
     * @param name Name of the array
     * @param platform_config Config parameter dictionary
     * @param column_names Columns to read
     * @param batch_size Read batch size
     * @param result_order Read result order
     * @return std::unique_ptr<SOMAArray> SOMAArray
     */
    static std::unique_ptr<SOMAArray> open(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::string_view name = "unnamed",
        std::map<std::string, std::string> platform_config = {},
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        std::string_view result_order = "auto",
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Open an array at the specified URI and return SOMAArray
     * object.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param ctx TileDB context
     * @param uri URI of the array
     * @param name Name of the array
     * @param column_names Columns to read
     * @param batch_size Read batch size
     * @param result_order Read result order
     * @return std::unique_ptr<SOMAArray> SOMAArray
     */
    static std::unique_ptr<SOMAArray> open(
        tiledb_query_type_t mode,
        std::shared_ptr<Context> ctx,
        std::string_view uri,
        std::string_view name = "unnamed",
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        std::string_view result_order = "auto",
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================
    /**
     * @brief Construct a new SOMAArray object
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param uri URI of the array
     * @param name name of the array
     * @param ctx TileDB context
     * @param column_names Columns to read
     * @param batch_size Batch size
     * @param result_order Result order
     * @param timestamp Timestamp
     */
    SOMAArray(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::string_view name,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        std::string_view result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    SOMAArray() = delete;
    SOMAArray(const SOMAArray&) = delete;
    SOMAArray(SOMAArray&&) = default;
    ~SOMAArray() = default;

    /**
     * Open the SOMAArray object.
     */
    void open(
        tiledb_query_type_t mode,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * Closes the SOMAArray object.
     */
    void close();

    /**
     * @brief Reset the state of this SOMAArray object to prepare for a
     * new query, while holding the array open.
     *
     * @param column_names
     * @param batch_size
     * @param result_order
     */
    void reset(
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        std::string_view result_order = "auto");

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
        mq_->select_point(dim, point);
    }

    /**
     * @brief Set the dimension slice using multiple points, with support for
     * partitioning.
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
        // Validate partition inputs
        if (partition_index >= partition_count) {
            throw TileDBSOMAError(fmt::format(
                "[SOMAArray] partition_index ({}) must be < "
                "partition_count "
                "({})",
                partition_index,
                partition_count));
        }

        if (partition_count > 1) {
            auto partition_size = points.size() / partition_count;
            auto start = partition_index * partition_size;

            // If this is the last partition, cover the rest of the points.
            if (partition_index == partition_count - 1) {
                partition_size = points.size() - start;
            }

            LOG_DEBUG(fmt::format(
                "[SOMAArray] set_dim_points partitioning: sizeof(T)={} "
                "dim={} "
                "index={} "
                "count={} "
                "range=[{}, {}] of {} points",
                sizeof(T),
                dim,
                partition_index,
                partition_count,
                start,
                start + partition_size - 1,
                points.size()));

            mq_->select_points(
                dim, tcb::span<T>{&points[start], partition_size});
        } else {
            mq_->select_points(dim, points);
        }
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
        LOG_DEBUG(
            fmt::format("[SOMAArray] set_dim_points: sizeof(T)={}", sizeof(T)));
        mq_->select_points(dim, points);
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
        mq_->select_ranges(dim, ranges);
    }

    /**
     * @brief Set a query condition.
     *
     * @param qc Query condition
     */
    void set_condition(QueryCondition& qc) {
        mq_->set_condition(qc);
    }

    /**
     * @brief Select columns names to query (dim and attr). If the
     * `if_not_empty` parameter is `true`, the column will be selected iff the
     * list of selected columns is empty. This prevents a `select_columns` call
     * from changing an empty list (all columns) to a subset of columns.
     *
     * @param names Vector of column names
     * @param if_not_empty Prevent changing an "empty" selection of all columns
     */
    void select_columns(
        const std::vector<std::string>& names, bool if_not_empty = false) {
        mq_->select_columns(names, if_not_empty);
    }

    /**
     * @brief Submit the query
     *
     */
    void submit();

    /**
     * @brief Read the next chunk of results from the query. If all results have
     * already been read, std::nullopt is returned.
     *
     * An example use model:
     *
     *   auto reader = SOMAArray::open(uri);
     *   reader->submit();
     *   while (auto batch = x_data->read_next()) {
     *       ...process batch ...
     *   }
     *
     * @return std::optional<std::shared_ptr<ArrayBuffers>>
     */
    std::optional<std::shared_ptr<ArrayBuffers>> read_next();

    /**
     * @brief Check if the query is complete.
     *
     * If `query_status_only` is true, return true if the query status is
     * complete.
     *
     * If `query_status_only` is false, return true if the query status
     * is complete or if the query is empty (no ranges have been added to the
     * query).
     *
     * @param query_status_only Query complete mode.
     * @return true if the query is complete, as described above
     */
    bool is_complete(bool query_status_only = false) {
        return mq_->is_complete(query_status_only);
    }

    /**
     * @brief Return true if `read_next` returned all results from the
     * query. The return value is false if the query was incomplete.
     *
     * @return True if last call to `read_next` returned all results of the
     * query
     */
    bool results_complete() {
        return mq_->results_complete();
    }

    /**
     * @brief Returns the total number of cells read so far, including any
     * previous incomplete queries.
     *
     * @return size_t Total number of cells read
     */
    size_t total_num_cells() {
        return mq_->total_num_cells();
    }

    /**
     * @brief Get the total number of unique cells in the array.
     *
     * @return uint64_t Total number of unique cells
     */
    uint64_t nnz();

    /**
     * @brief Get the schema of the array.
     *
     * @return std::shared_ptr<ArraySchema> Schema
     */
    std::shared_ptr<ArraySchema> schema() {
        return mq_->schema();
    }

    /**
     * @brief Get the capacity of each dimension.
     *
     * @return A vector with length equal to the number of dimensions; each
     * value in the vector is the capcity of each dimension.
     */
    std::vector<int64_t> shape();

    /**
     * Set metadata key-value items to an open array. The array must
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
        const void* value);

    /**
     * Deletes a metadata key-value item from an open array. The array must
     * be opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be deleted.
     *
     * @note The writes will take effect only upon closing the array.
     *
     * @note If the key does not exist, this will take no effect
     *     (i.e., the function will not error out).
     */
    void delete_metadata(const std::string& key);

    /**
     * @brief Given a key, retrieve the associated value datatype, number of
     * values, and value in binary form. The array must be opened in READ mode,
     * otherwise the function will error out.
     *
     * The value may consist of more than one items of the same datatype. Keys
     * that do not exist in the metadata will be return NULL for the value.
     *
     * @param key The key of the metadata item to be retrieved. UTF-8 encodings
     *     are acceptable.
     * @return MetadataValue (std::tuple<std::string, tiledb_datatype_t,
     * uint32_t, const void*>)
     */
    MetadataValue get_metadata(const std::string& key) const;

    /**
     * @brief Given an index, retrieve the associated value datatype, number of
     * values, and value in binary form. The array must be opened in READ mode,
     * otherwise the function will error out.
     *
     * @param index The index used to get the metadata.
     * @return MetadataValue (std::tuple<std::string, tiledb_datatype_t,
     * uint32_t, const void*>)
     */
    MetadataValue get_metadata(uint64_t index) const;

    /**
     * Checks if key exists in metadata from an open array. The array must
     * be opened in READ mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be checked. UTF-8 encodings
     *     are acceptable.
     * @return true if the key exists, else false.
     */
    bool has_metadata(const std::string& key);

    /**
     * Returns then number of metadata items in an open array. The array must
     * be opened in READ mode, otherwise the function will error out.
     */
    uint64_t metadata_num() const;

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMAArray URI
    std::string uri_;

    // Batch size
    std::string batch_size_;

    // Result order
    std::string result_order_;

    // Read timestamp range (start, end)
    std::optional<std::pair<uint64_t, uint64_t>> timestamp_;

    // Managed query for the array
    std::unique_ptr<ManagedQuery> mq_;

    // Array associated with mq_
    std::shared_ptr<Array> arr_;

    // True if this is the first call to read_next()
    bool first_read_next_ = true;

    // True if the query was submitted
    bool submitted_ = false;

    // Unoptimized method for computing nnz() (issue `count_cells` query)
    uint64_t nnz_slow();
};

}  // namespace tiledbsoma

#endif  // SOMA_ARRAY
