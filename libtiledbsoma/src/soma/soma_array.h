/**
 * @file   soma_array.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022-2024 TileDB, Inc.
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
#include "../utils/arrow_adapter.h"
#include "enums.h"
#include "logger_public.h"
#include "managed_query.h"
#include "soma_object.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAArray : public SOMAObject {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAArray object at the given URI.
     *
     * @param ctx SOMAContext
     * @param uri URI to create the SOMAArray
     * @param schema TileDB ArraySchema
     * @param soma_type SOMADataFrame, SOMADenseNDArray, or
     * SOMASparseNDArray
     */
    static std::unique_ptr<SOMAArray> create(
        std::shared_ptr<SOMAContext> ctx,
        std::string_view uri,
        ArraySchema schema,
        std::string soma_type,
        std::optional<TimestampRange> timestamp = std::nullopt);

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
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open an array at the specified URI and return SOMAArray
     * object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx SOMAContext
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
        std::shared_ptr<SOMAContext> ctx,
        std::string_view name = "unnamed",
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<TimestampRange> timestamp = std::nullopt);

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
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Construct a new SOMAArray object
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx SOMAContext
     * @param name name of the array
     * @param column_names Columns to read
     * @param batch_size Batch size
     * @param result_order Result order
     * @param timestamp Timestamp
     */
    SOMAArray(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::string_view name,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        ResultOrder result_order,
        std::optional<TimestampRange> timestamp = std::nullopt);

    SOMAArray(const SOMAArray& other)
        : uri_(other.uri_)
        , name_(other.name_)
        , ctx_(other.ctx_)
        , batch_size_(other.batch_size_)
        , result_order_(other.result_order_)
        , metadata_(other.metadata_)
        , timestamp_(other.timestamp_)
        , mq_(std::make_unique<ManagedQuery>(
              other.arr_, other.ctx_->tiledb_ctx(), other.name_))
        , arr_(other.arr_)
        , meta_cache_arr_(other.meta_cache_arr_)
        , first_read_next_(other.first_read_next_)
        , submitted_(other.submitted_)
        , array_buffer_(other.array_buffer_) {
        fill_metadata_cache();
    }

    SOMAArray(
        std::shared_ptr<SOMAContext> ctx,
        std::shared_ptr<Array> arr,
        std::optional<TimestampRange> timestamp);

    SOMAArray(SOMAArray&&) = default;

    SOMAArray(const SOMAObject& other)
        : SOMAObject(other) {
    }

    SOMAArray() = delete;
    ~SOMAArray() = default;

    /**
     * @brief Get URI of the SOMAArray.
     *
     * @return std::string URI
     */
    const std::string uri() const;

    /**
     * @brief Get context of the SOMAArray.
     *
     * @return SOMAContext
     */
    std::shared_ptr<SOMAContext> ctx();

    /**
     * Open the SOMAArray object.
     *
     * @param mode read or write
     * @param timestamp Timestamp
     */
    void open(
        OpenMode mode, std::optional<TimestampRange> timestamp = std::nullopt);

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
     * Get whether the SOMAArray was open in read or write mode.
     *
     * @return OpenMode
     */
    OpenMode mode() const {
        return mq_->query_type() == TILEDB_READ ? OpenMode::read :
                                                  OpenMode::write;
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

    Enumeration extend_enumeration(
        ArrowSchema* value_schema,
        ArrowArray* value_array,
        ArrowSchema* index_schema,
        ArrowArray* index_array);

    /**
     * @brief Set the write buffers for a single column.
     *
     * @param name Name of the column
     * @param num_elems Number of elements to write
     * @param data Pointer to the beginning of the data buffer
     * @param offsets Optional pointer to the beginning of the offsets buffer
     * @param validity Optional pointer to the beginning of the validities
     * buffer
     */
    void set_column_data(
        std::string_view name,
        uint64_t num_elems,
        const void* data,
        uint64_t* offsets = nullptr,
        uint8_t* validity = nullptr);

    /**
     * @brief Set the write buffers for string or binary with 32-bit offsets
     * (as opposed to large string or large binary with 64-bit offsets).
     *
     * @param name Name of the column
     * @param num_elems Number of elements to write
     * @param data Pointer to the beginning of the data buffer
     * @param offsets Pointer to the beginning of the offsets buffer
     * @param validity Optional pointer to the beginning of the validities
     * buffer
     */
    void set_column_data(
        std::string_view name,
        uint64_t num_elems,
        const void* data,
        uint32_t* offsets,
        uint8_t* validity = nullptr);

    /**
     * @brief Set the write buffers for an Arrow Table or Batch as represented
     * by an ArrowSchema and ArrowArray.
     *
     * @param arrow_schema
     * @param arrow_array
     */
    void set_array_data(
        std::unique_ptr<ArrowSchema> arrow_schema,
        std::unique_ptr<ArrowArray> arrow_array);

    /**
     * @brief Write ArrayBuffers data to the array after setting write buffers.
     *
     * An example use model:
     *
     *   auto array = SOMAArray::open(TILEDB_WRITE, uri);
     *   array.set_array_data(
     *      std::make_unique<ArrowSchema>(arrow_schema),
     *      std::make_unique<ArrowArray>(arrow_array));
     *   array.write();
     *   array.close();
     */
    void write(bool sort_coords = true);

    /**
     * @brief Consolidates and vacuums fragment metadata and commit files.
     *
     * @param modes List of modes to apply. By default, apply to fragment_meta
     * and commits
     */
    void consolidate_and_vacuum(
        std::vector<std::string> modes = {"fragment_meta", "commits"});

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
     * @brief Get the TileDB ArraySchema. This should eventually
     * be removed in lieu of arrow_schema below.
     *
     * @return std::shared_ptr<ArraySchema> Schema
     */
    std::shared_ptr<ArraySchema> tiledb_schema() const {
        return mq_->schema();
    }

    /**
     * @brief Get the Arrow schema of the array.
     *
     * @return std::unique_ptr<ArrowSchema> Schema
     */
    std::unique_ptr<ArrowSchema> arrow_schema() const {
        return ArrowAdapter::arrow_schema_from_tiledb_array(
            ctx_->tiledb_ctx(), arr_);
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
     * Retrieves the non-empty domain from the array. This is the union of the
     * non-empty domains of the array fragments.
     */
    template <typename T>
    std::pair<T, T> non_empty_domain(const std::string& name) {
        try {
            return arr_->non_empty_domain<T>(name);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        }
    }

    /**
     * Retrieves the non-empty domain from the array on the given dimension.
     * This is the union of the non-empty domains of the array fragments.
     * Applicable only to var-sized dimensions.
     */
    std::pair<std::string, std::string> non_empty_domain_var(
        const std::string& name) {
        try {
            return arr_->non_empty_domain_var(name);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        }
    }

    /**
     * Returns the domain of the given dimension.
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    template <typename T>
    std::pair<T, T> domain(const std::string& name) const {
        return arr_->schema().domain().dimension(name).domain<T>();
    }

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
        std::optional<TimestampRange> timestamp);

    /**
     * Return optional timestamp pair SOMAArray was opened with.
     */
    std::optional<TimestampRange> timestamp();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    uint64_t _get_max_capacity(tiledb_datatype_t index_type);

    template <typename ValueType>
    Enumeration _extend_and_evolve_schema(
        ArrowArray* value_array,
        ArrowSchema* index_schema,
        ArrowArray* index_array) {
        std::string column_name = index_schema->name;
        auto disk_index_type = tiledb_schema()->attribute(column_name).type();
        auto enmr = ArrayExperimental::get_enumeration(
            *ctx_->tiledb_ctx(), *arr_, column_name);
        uint64_t max_capacity = SOMAArray::_get_max_capacity(disk_index_type);

        const void* data;
        uint64_t num_elems = value_array->length;
        if (value_array->n_buffers == 3) {
            data = value_array->buffers[2];
        } else {
            data = value_array->buffers[1];
        }

        std::vector<ValueType> enums_in_write(
            (ValueType*)data, (ValueType*)data + num_elems);
        auto enums_existing = enmr.as_vector<ValueType>();
        std::vector<ValueType> extend_values;
        for (auto enum_val : enums_in_write) {
            if (std::find(
                    enums_existing.begin(), enums_existing.end(), enum_val) ==
                enums_existing.end()) {
                extend_values.push_back(enum_val);
            }
        }

        if (extend_values.size() != 0) {
            auto free_capacity = max_capacity - enums_existing.size();
            if (free_capacity < extend_values.size()) {
                throw TileDBSOMAError(
                    "Cannot extend enumeration; reached maximum capacity");
            }
            ArraySchemaEvolution se(*ctx_->tiledb_ctx());
            auto extended_enmr = enmr.extend(extend_values);
            se.extend_enumeration(extended_enmr);
            se.array_evolve(uri_);
            SOMAArray::_remap_indexes(
                column_name,
                extended_enmr,
                enums_in_write,
                index_schema,
                index_array);

            index_schema->format = ArrowAdapter::to_arrow_format(
                                       disk_index_type)
                                       .data();

            return extended_enmr;
        }

        return enmr;
    }

    Enumeration _extend_and_evolve_schema_str(
        ArrowSchema* value_schema,
        ArrowArray* value_array,
        ArrowSchema* index_schema,
        ArrowArray* index_array);

    template <typename UserType>
    void _cast_column(
        ArrowSchema* orig_column_schema,
        ArrowArray* orig_column_array,
        // ArrowSchema* new_column_schema,
        ArrowArray* new_column_array) {
        tiledb_datatype_t disk_type;
        std::string name(orig_column_schema->name);
        if (tiledb_schema()->has_attribute(name)) {
            disk_type = tiledb_schema()->attribute(name).type();
        } else {
            disk_type = tiledb_schema()->domain().dimension(name).type();
        }

        std::cout << " " << tiledb::impl::type_to_str(disk_type) << std::endl;

        switch (disk_type) {
            case TILEDB_STRING_ASCII:
            case TILEDB_STRING_UTF8:
            case TILEDB_CHAR:
                break;
            case TILEDB_BOOL:
            case TILEDB_INT8:
                return SOMAArray::_cast_column_aux<UserType, int8_t>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_UINT8:
                return SOMAArray::_cast_column_aux<UserType, uint8_t>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_INT16:
                return SOMAArray::_cast_column_aux<UserType, int16_t>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_UINT16:
                return SOMAArray::_cast_column_aux<UserType, uint16_t>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_INT32:
                return SOMAArray::_cast_column_aux<UserType, int32_t>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_UINT32:
                return SOMAArray::_cast_column_aux<UserType, uint32_t>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_INT64:
            case TILEDB_DATETIME_YEAR:
            case TILEDB_DATETIME_MONTH:
            case TILEDB_DATETIME_WEEK:
            case TILEDB_DATETIME_DAY:
            case TILEDB_DATETIME_HR:
            case TILEDB_DATETIME_MIN:
            case TILEDB_DATETIME_SEC:
            case TILEDB_DATETIME_MS:
            case TILEDB_DATETIME_US:
            case TILEDB_DATETIME_NS:
            case TILEDB_DATETIME_PS:
            case TILEDB_DATETIME_FS:
            case TILEDB_DATETIME_AS:
            case TILEDB_TIME_HR:
            case TILEDB_TIME_MIN:
            case TILEDB_TIME_SEC:
            case TILEDB_TIME_MS:
            case TILEDB_TIME_US:
            case TILEDB_TIME_NS:
            case TILEDB_TIME_PS:
            case TILEDB_TIME_FS:
            case TILEDB_TIME_AS:
                return SOMAArray::_cast_column_aux<UserType, int64_t>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_UINT64:
                return SOMAArray::_cast_column_aux<UserType, uint64_t>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_FLOAT32:
                return SOMAArray::_cast_column_aux<UserType, float>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            case TILEDB_FLOAT64:
                return SOMAArray::_cast_column_aux<UserType, double>(
                    // orig_column_schema,
                    orig_column_array,
                    // new_column_schema,
                    new_column_array);
            default:
                throw TileDBSOMAError(
                    "Saw invalid TileDB disk type when attempting to cast "
                    "table: " +
                    tiledb::impl::type_to_str(disk_type));
        }
    }

    template <typename UserType, typename DiskType>
    void _cast_column_aux(
        // ArrowSchema* orig_column_schema,
        ArrowArray* orig_column_array,
        // ArrowSchema* new_column_schema,
        ArrowArray* new_column_array) {
        UserType* buf;
        if (orig_column_array->n_buffers == 3) {
            buf = (UserType*)orig_column_array->buffers[2] +
                  orig_column_array->offset;
        } else {
            buf = (UserType*)orig_column_array->buffers[1] +
                  orig_column_array->offset;
        }
        std::vector<UserType> original_values(
            buf, buf + orig_column_array->length);
        std::vector<DiskType> casted_values;
        for (auto val : original_values) {
            casted_values.push_back(val);
        }

        if (orig_column_array->n_buffers == 3) {
            new_column_array->buffers[2] = malloc(
                sizeof(DiskType) * casted_values.size());
            std::memcpy(
                (void*)new_column_array->buffers[2],
                casted_values.data(),
                sizeof(DiskType) * casted_values.size());
        } else {
            new_column_array->buffers[1] = malloc(
                sizeof(DiskType) * casted_values.size());
            std::memcpy(
                (void*)new_column_array->buffers[1],
                casted_values.data(),
                sizeof(DiskType) * casted_values.size());
        }
    }

    template <typename ValueType>
    void _remap_indexes(
        std::string column_name,
        Enumeration extended_enmr,
        std::vector<ValueType> enums_in_write,
        ArrowSchema* index_schema,
        ArrowArray* index_array) {
        auto user_index_type = ArrowAdapter::to_tiledb_format(
            index_schema->format);
        switch (user_index_type) {
            case TILEDB_INT8:
                SOMAArray::_remap_indexes_aux<ValueType, int8_t>(
                    column_name, extended_enmr, enums_in_write, index_array);
                break;
            case TILEDB_UINT8:
                SOMAArray::_remap_indexes_aux<ValueType, uint8_t>(
                    column_name, extended_enmr, enums_in_write, index_array);
                break;
            case TILEDB_INT16:
                SOMAArray::_remap_indexes_aux<ValueType, int16_t>(
                    column_name, extended_enmr, enums_in_write, index_array);
                break;
            case TILEDB_UINT16:
                SOMAArray::_remap_indexes_aux<ValueType, uint16_t>(
                    column_name, extended_enmr, enums_in_write, index_array);
                break;
            case TILEDB_INT32:
                SOMAArray::_remap_indexes_aux<ValueType, int32_t>(
                    column_name, extended_enmr, enums_in_write, index_array);
                break;
            case TILEDB_UINT32:
                SOMAArray::_remap_indexes_aux<ValueType, uint32_t>(
                    column_name, extended_enmr, enums_in_write, index_array);
                break;
            case TILEDB_INT64:
                SOMAArray::_remap_indexes_aux<ValueType, int64_t>(
                    column_name, extended_enmr, enums_in_write, index_array);
                break;
            case TILEDB_UINT64:
                SOMAArray::_remap_indexes_aux<ValueType, uint64_t>(
                    column_name, extended_enmr, enums_in_write, index_array);
                break;
            default:
                throw TileDBSOMAError(
                    "Saw invalid enumeration index type when trying to extend"
                    "enumeration");
        }
    }

    template <typename ValueType, typename IndexType>
    void _remap_indexes_aux(
        std::string column_name,
        Enumeration extended_enmr,
        std::vector<ValueType> enums_in_write,
        ArrowArray* index_array) {
        IndexType* idxbuf;
        if (index_array->n_buffers == 3) {
            idxbuf = (IndexType*)index_array->buffers[2];
        } else {
            idxbuf = (IndexType*)index_array->buffers[1];
        }

        auto enmr_vec = extended_enmr.as_vector<ValueType>();
        auto beg = enmr_vec.begin();
        auto end = enmr_vec.end();
        std::vector<IndexType> original_indexes(
            idxbuf, idxbuf + index_array->length);
        std::vector<IndexType> shifted_indexes;
        for (auto i : original_indexes) {
            auto it = std::find(beg, end, enums_in_write[i]);
            shifted_indexes.push_back(it - beg);
        }

        auto attr = tiledb_schema()->attribute(column_name);
        switch (attr.type()) {
            case TILEDB_INT8:
                return SOMAArray::_cast_shifted_indexes<IndexType, int8_t>(
                    shifted_indexes, index_array);
            case TILEDB_UINT8:
                return SOMAArray::_cast_shifted_indexes<IndexType, uint8_t>(
                    shifted_indexes, index_array);
            case TILEDB_INT16:
                return SOMAArray::_cast_shifted_indexes<IndexType, int16_t>(
                    shifted_indexes, index_array);
            case TILEDB_UINT16:
                return SOMAArray::_cast_shifted_indexes<IndexType, uint16_t>(
                    shifted_indexes, index_array);
            case TILEDB_INT32:
                return SOMAArray::_cast_shifted_indexes<IndexType, int32_t>(
                    shifted_indexes, index_array);
            case TILEDB_UINT32:
                return SOMAArray::_cast_shifted_indexes<IndexType, uint32_t>(
                    shifted_indexes, index_array);
            case TILEDB_INT64:
                return SOMAArray::_cast_shifted_indexes<IndexType, int64_t>(
                    shifted_indexes, index_array);
            case TILEDB_UINT64:
                return SOMAArray::_cast_shifted_indexes<IndexType, uint64_t>(
                    shifted_indexes, index_array);
            default:
                throw TileDBSOMAError(
                    "Saw invalid enumeration index type when trying to extend"
                    "enumeration");
        }
    }

    template <typename UserIndexType, typename DiskIndexType>
    void _cast_shifted_indexes(
        std::vector<UserIndexType> shifted_indexes, ArrowArray* index_array) {
        std::vector<DiskIndexType> casted_indexes;
        for (auto i : shifted_indexes) {
            casted_indexes.push_back(i);
        }

        if (index_array->n_buffers == 3) {
            index_array->buffers[2] = malloc(
                sizeof(DiskIndexType) * casted_indexes.size());
            std::memcpy(
                (void*)index_array->buffers[2],
                casted_indexes.data(),
                sizeof(DiskIndexType) * casted_indexes.size());
        } else {
            index_array->buffers[1] = malloc(
                sizeof(DiskIndexType) * casted_indexes.size());
            std::memcpy(
                (void*)index_array->buffers[1],
                casted_indexes.data(),
                sizeof(DiskIndexType) * casted_indexes.size());
        }
    }

    // Helper function to cast Boolean of bits (Arrow) to uint8 (TileDB)
    void _cast_bit_to_uint8(ArrowSchema* arrow_schema, ArrowArray* arrow_array);

    // Helper function for set_column_data
    std::shared_ptr<ColumnBuffer> _setup_column_data(std::string_view name);

    // Fills the metadata cache upon opening the array.
    void fill_metadata_cache();

    // Helper function for set_array_data
    ArrowTable _cast_table(
        std::unique_ptr<ArrowSchema> arrow_schema,
        std::unique_ptr<ArrowArray> arrow_array);

    // SOMAArray URI
    std::string uri_;

    // SOMAArray name for debugging
    std::string_view name_;

    // SOMA context
    std::shared_ptr<SOMAContext> ctx_;

    // Batch size
    std::string batch_size_;

    // Result order
    ResultOrder result_order_;

    // Metadata cache
    std::map<std::string, MetadataValue> metadata_;

    // Read timestamp range (start, end)
    std::optional<TimestampRange> timestamp_;

    // Managed query for the array
    std::unique_ptr<ManagedQuery> mq_;

    // Array associated with mq_
    std::shared_ptr<Array> arr_;

    // Array associated with metadata_. Metadata values need to be accessible in
    // write mode as well. We need to keep this read-mode array alive in order
    // for the metadata value pointers in the cache to be accessible
    std::shared_ptr<Array> meta_cache_arr_;

    // True if this is the first call to read_next()
    bool first_read_next_ = true;

    // True if the query was submitted
    bool submitted_ = false;

    // Unoptimized method for computing nnz() (issue `count_cells` query)
    uint64_t nnz_slow();

    // ArrayBuffers to hold ColumnBuffers alive when submitting to write query
    std::shared_ptr<ArrayBuffers> array_buffer_ = nullptr;
};

}  // namespace tiledbsoma

#endif  // SOMA_ARRAY
