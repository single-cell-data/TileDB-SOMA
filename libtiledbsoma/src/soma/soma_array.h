/**
 * @file   soma_array.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022-2023 TileDB, Inc.
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
 *   This declares the SOMAArray class.
 */

#ifndef SOMA_ARRAY
#define SOMA_ARRAY

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <future>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>
#include "enums.h"
#include "logger_public.h"
#include "managed_query.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAArray {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAArray object at the given URI.
     *
     * @param ctx TileDB context
     * @param uri URI to create the SOMAArray
     * @param schema TileDB ArraySchema
     * @param soma_type SOMADataFrame, SOMADenseNDArray, or
     * SOMASparseNDArray
     */
    static void create(
        std::shared_ptr<Context> ctx,
        std::string_view uri,
        ArraySchema schema,
        std::string soma_type);

    /**
     * @brief Open an array at the specified URI and return SOMAArray
     * object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param name Name of the array
     * @param platform_config Config parameter dictionary
     * @param column_names Columns to read
     * @param batch_size Read batch size
     * @param result_order Read result order: automatic (default), rowmajor,
     * or colmajor
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::unique_ptr<SOMAArray> SOMAArray
     */
    static std::unique_ptr<SOMAArray> open(
        OpenMode mode,
        std::string_view uri,
        std::string_view name = "unnamed",
        std::map<std::string, std::string> platform_config = {},
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Open an array at the specified URI and return SOMAArray
     * object.
     *
     * @param mode read or write
     * @param ctx TileDB context
     * @param uri URI of the array
     * @param name Name of the array
     * @param column_names Columns to read
     * @param batch_size Read batch size
     * @param result_order Read result order: automatic (default), rowmajor,
     * or colmajor
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::unique_ptr<SOMAArray> SOMAArray
     */
    static std::unique_ptr<SOMAArray> open(
        OpenMode mode,
        std::shared_ptr<Context> ctx,
        std::string_view uri,
        std::string_view name = "unnamed",
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMAArray object
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param name name of the array
     * @param platform_config Config parameter dictionary
     * @param column_names Columns to read
     * @param batch_size Batch size
     * @param result_order Result order
     * @param timestamp Timestamp
     */
    SOMAArray(
        OpenMode mode,
        std::string_view uri,
        std::string_view name,
        std::map<std::string, std::string> platform_config,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        ResultOrder result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Construct a new SOMAArray object
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param name name of the array
     * @param ctx TileDB context
     * @param column_names Columns to read
     * @param batch_size Batch size
     * @param result_order Read result order: automatic (default), rowmajor,
     * or colmajor
     * @param timestamp Timestamp
     */
    SOMAArray(
        OpenMode mode,
        std::string_view uri,
        std::string_view name,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        ResultOrder result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    SOMAArray() = delete;
    SOMAArray(const SOMAArray&) = delete;
    SOMAArray(SOMAArray&&) = default;
    ~SOMAArray() = default;

    /**
     * @brief Get URI of the SOMAArray.
     *
     * @return std::string URI
     */
    const std::string& uri() const;

    /**
     * @brief Get Ctx of the SOMAArray.
     *
     * @return std::shared_ptr<Context>
     */
    std::shared_ptr<Context> ctx();

    /**
     * Open the SOMAArray object.
     *
     * @param mode read or write
     * @param timestamp Timestamp
     */
    void open(
        OpenMode mode,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * Close the SOMAArray object.
     */
    void close();

    /**
     * Check if the SOMAArray is open.
     *
     * @return bool true if open
     */
    bool is_open() const {
        return arr_->is_open();
    }

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
        ResultOrder result_order = ResultOrder::automatic);

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
        // Validate partition inputs
        if (partition_index >= partition_count) {
            // TODO this use to be formatted with fmt::format which is part of
            // internal header spd/log/fmt/fmt.h and should not be used.
            // In C++20, this can be replaced with std::format.
            std::ostringstream err;
            err << "[SOMAArray] partition_index (" << partition_index
                << ") must be < partition_count (" << partition_count;
            throw TileDBSOMAError(err.str());
        }

        if (partition_count > 1) {
            auto partition_size = points.size() / partition_count;
            auto start = partition_index * partition_size;

            // If this is the last partition, cover the rest of the points.
            if (partition_index == partition_count - 1) {
                partition_size = points.size() - start;
            }

            // TODO this use to be formatted with fmt::format which is part of
            // internal header spd/log/fmt/fmt.h and should not be used.
            // In C++20, this can be replaced with std::format.
            std::ostringstream log_dbg;
            log_dbg << "[SOMAArray] set_dim_points partitioning:"
                    << " sizeof(T)=" << sizeof(T) << " dim=" << dim
                    << " index=" << partition_index
                    << " count=" << partition_count << " range =[" << start
                    << ", " << start + partition_size - 1 << "] of "
                    << points.size() << "points";
            LOG_DEBUG(log_dbg.str());

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
            "[SOMAArray] set_dim_points: sizeof(T)=" +
            std::to_string(sizeof(T)));
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
     * `if_not_empty` parameter is `true`, the column will be selected iff
     * the list of selected columns is empty. This prevents a
     * `select_columns` call from changing an empty list (all columns) to a
     * subset of columns.
     *
     * @param names Vector of column names
     * @param if_not_empty Prevent changing an "empty" selection of all
     * columns
     */
    void select_columns(
        const std::vector<std::string>& names, bool if_not_empty = false) {
        mq_->select_columns(names, if_not_empty);
    }

    /**
     * @brief Returns the column names set by the query.
     *
     * @return std::vector<std::string>
     */
    std::vector<std::string> column_names() {
        return mq_->column_names();
    }

    /**
     * @brief Returns the result order set by the query.
     *
     * @return ResultOrder
     */
    ResultOrder result_order() {
        return result_order_;
    }

    /**
     * @brief Read the next chunk of results from the query. If all results
     * have already been read, std::nullopt is returned.
     *
     * An example use model:
     *
     *   auto reader = SOMAArray::open(TILEDB_READ, uri);
     *   while (auto batch = x_data->read_next()) {
     *       ...process batch ...
     *   }
     *
     * @return std::optional<std::shared_ptr<ArrayBuffers>>
     */
    std::optional<std::shared_ptr<ArrayBuffers>> read_next();

    /**
     * @brief Set the write data for a column.
     *
     * @param column_name Column name
     * @param buff Buffer array pointer with elements of the column type.
     * @param nelements Number of array elements in buffer
     */
    void set_column_data(
        std::string_view column_name,
        std::shared_ptr<ColumnBuffer> column_buffer) {
        mq_->set_column_data(std::string(column_name), column_buffer);
    }

    /**
     * @brief Write ArrayBuffers data to the array.
     *
     * An example use model:
     *
     *   auto writer = SOMAArray::open(TILEDB_WRITE, uri);
     *
     *   std::vector<int> att{0, 1, 2, 3, 4, 5};
     *   std::vector<int> dim{0, 1, 2, 3, 4, 5};
     *
     *   auto schema = *soma_array->schema();
     *   auto array_buffer = std::make_shared<ArrayBuffers>();
     *   array_buffer->emplace("att", ColumnBuffer::create(schema, "att",
     * att)); array_buffer->emplace("dim", ColumnBuffer::create(schema,
     * "dim", dim));
     *
     *   std::vector<int> x(10, 1);
     *   writer->write(array_buffer);
     *   writer->close();
     *
     * @param buffers The ArrayBuffers to write to the array
     */
    void write(std::shared_ptr<ArrayBuffers> buffers);

    /**
     * @brief Check if the query is complete.
     *
     * If `query_status_only` is true, return true if the query status is
     * complete.
     *
     * If `query_status_only` is false, return true if the query status
     * is complete or if the query is empty (no ranges have been added to
     * the query).
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
     * @brief Return the total number of cells read so far, including any
     * previous incomplete queries.
     *
     * @return size_t Total number of cells read
     */
    size_t total_num_cells() {
        return mq_->total_num_cells();
    }

    /**
     * @brief Return whether next read is the initial read, or a subsequent
     * read of a previous incomplete query.
     *
     * @return bool Logical value if initial read or not
     */
    bool is_initial_read() {
        return first_read_next_;
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
    std::shared_ptr<ArraySchema> schema() const {
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
     * @brief Get the number of dimensions.
     *
     * @return uint64_t Number of dimensions.
     */
    uint64_t ndim() const;

    /**
     * @brief Get the name of each dimensions.
     *
     * @return std::vector<std::string> Name of each dimensions.
     */
    std::vector<std::string> dimension_names() const;

    /**
     * @brief Get the mapping of attributes to Enumerations.
     *
     * @return std::map<std::string, Enumeration>
     */
    std::map<std::string, Enumeration> get_attr_to_enum_mapping();

    /**
     * @brief Get the Enumeration name associated with the given Attr.
     *
     * @return std::optional<std::string> The enumeration name if one exists.
     */
    std::optional<std::string> get_enum_label_on_attr(std::string attr_name);

    /**
     * @brief Check if the given attribute has an associated enumeration.
     *
     * @return bool
     */
    bool attr_has_enum(std::string attr_name);

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
     * Delete a metadata key-value item from an open array. The array must
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
     * @brief Given a key, get the associated value datatype, number of
     * values, and value in binary form. The array must be opened in READ
     mode,
     * otherwise the function will error out.
     *
     * The value may consist of more than one items of the same datatype.
     Keys
     * that do not exist in the metadata will be return NULL for the value.
     *
     * **Example:**
     * @code{.cpp}
     * // Open the array for reading
     * tiledbsoma::SOMAArray soma_array = SOMAArray::open(TILEDB_READ,
     "s3://bucket-name/group-name");
     * tiledbsoma::MetadataValue meta_val = soma_array->get_metadata("key");
     * std::string key = std::get<MetadataInfo::key>(meta_val);
     * tiledb_datatype_t dtype = std::get<MetadataInfo::dtype>(meta_val);
     * uint32_t num = std::get<MetadataInfo::num>(meta_val);
     * const void* value = *((const
     int32_t*)std::get<MetadataInfo::value>(meta_val));
     * @endcode
     *
     * @param key The key of the metadata item to be retrieved. UTF-8
     encodings
     *     are acceptable.
     * @return MetadataValue (std::tuple<std::string, tiledb_datatype_t,
     * uint32_t, const void*>)
     */
    std::optional<MetadataValue> get_metadata(const std::string& key);

    /**
     * Get a mapping of all metadata keys with its associated value datatype,
     * number of values, and value in binary form.
     *
     * @return std::map<std::string, MetadataValue>
     */
    std::map<std::string, MetadataValue> get_metadata();

    /**
     * Check if the key exists in metadata from an open array. The array
     * must be opened in READ mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be checked. UTF-8
     * encodings are acceptable.
     * @return true if the key exists, else false.
     */
    bool has_metadata(const std::string& key);

    /**
     * Return then number of metadata items in an open array. The array must
     * be opened in READ mode, otherwise the function will error out.
     */
    uint64_t metadata_num() const;

    /**
     * Validates input parameters before opening array.
     */
    void validate(
        OpenMode mode,
        std::string_view name,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp);

    /**
     * Return optional timestamp pair SOMAArray was opened with.
     */
    std::optional<std::pair<uint64_t, uint64_t>> timestamp();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    /**
     * Fills the metadata cache upon opening the array.
     */
    void fill_metadata_cache();

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMAArray URI
    std::string uri_;

    // Batch size
    std::string batch_size_;

    // Result order
    ResultOrder result_order_;

    // Metadata cache
    std::map<std::string, MetadataValue> metadata_;

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
