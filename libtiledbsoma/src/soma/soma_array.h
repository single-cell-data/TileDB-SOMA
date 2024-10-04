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

// This enables some code deduplication between core domain, core current
// domain, and core non-empty domain.
enum class Domainish {
    kind_core_domain = 0,
    kind_core_current_domain = 1,
    kind_non_empty_domain = 2
};

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
        , submitted_(other.submitted_) {
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
     * Return a new SOMAArray with the given mode at the current Unix timestamp.
     *
     * @param mode if the OpenMode is not given, If the SOMAObject was opened in
     * READ mode, reopen it in WRITE mode and vice versa
     * @param timestamp Timestamp
     */
    std::unique_ptr<SOMAArray> reopen(
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
     * @brief Get the number of dimensions.
     *
     * @return uint64_t Number of dimensions.
     */
    uint64_t ndim() const;

    /**
     * @brief Get the name of each dimension.
     *
     * @return std::vector<std::string> Name of each dimension.
     */
    std::vector<std::string> dimension_names() const;

    /**
     * @brief Sees if the array has a dimension of the given name.
     *
     * @return bool
     */
    bool has_dimension_name(const std::string& name) const;

    /**
     * @brief Get the name of each attribute.
     *
     * @return std::vector<std::string> Name of each attribute.
     */
    std::vector<std::string> attribute_names() const;

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
     * @param force A boolean toggle to suppress internal checks, defaults to
     *     false.
     *
     * @note The writes will take effect only upon closing the array.
     */
    void set_metadata(
        const std::string& key,
        tiledb_datatype_t value_type,
        uint32_t value_num,
        const void* value,
        bool force = false);

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
     * mode, otherwise the function will error out.
     *
     * The value may consist of more than one items of the same datatype.
     * Keys that do not exist in the metadata will be return NULL for the value.
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
     * encodings are acceptable.
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

    /**
     * Retrieves the non-empty domain from the array. This is the union of the
     * non-empty domains of the array fragments.
     */
    template <typename T>
    std::pair<T, T> non_empty_domain_slot(const std::string& name) const {
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
    std::pair<std::string, std::string> non_empty_domain_slot_var(
        const std::string& name) const {
        try {
            return arr_->non_empty_domain_var(name);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        }
    }

    /**
     * Exposed for testing purposes within this library.
     * Not for use by Python/R.
     */
    CurrentDomain get_current_domain_for_test() const {
        return _get_current_domain();
    }

    /**
     * @brief Returns true if the array has a non-empty current domain, else
     * false.  Note that at the core level it's "current domain" for all arrays;
     * at the SOMA-API level it's "upgraded_shape" for SOMASparseNDArray and
     * SOMADenseNDArray, and "upgraded_domain" for SOMADataFrame; here
     * we use the core language and defer to Python/R to conform to
     * SOMA-API syntax.
     */
    bool has_current_domain() const {
        return !_get_current_domain().is_empty();
    }

    /**
     * Returns the core current domain at the given dimension.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     *   - soma maxdomain is core domain
     *   - core current domain is not accessed at the soma level
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    template <typename T>
    std::pair<T, T> _core_current_domain_slot(const std::string& name) const {
        if (std::is_same_v<T, std::string>) {
            throw std::runtime_error(
                "SOMAArray::soma_domain_slot: template-specialization "
                "failure.");
        }
        CurrentDomain current_domain = _get_current_domain();
        if (current_domain.is_empty()) {
            throw TileDBSOMAError(
                "_core_current_domain_slot: internal coding error");
        }
        if (current_domain.type() != TILEDB_NDRECTANGLE) {
            throw TileDBSOMAError(
                "_core_current_domain_slot: found non-rectangle type");
        }
        NDRectangle ndrect = current_domain.ndrectangle();

        // Convert from two-element array (core API) to pair (tiledbsoma API)
        std::array<T, 2> arr = ndrect.range<T>(name);
        return std::pair<T, T>(arr[0], arr[1]);
    }

    std::pair<std::string, std::string> _core_current_domain_slot_string(
        const std::string& name) const {
        CurrentDomain current_domain = _get_current_domain();
        if (current_domain.is_empty()) {
            throw TileDBSOMAError(
                "_core_current_domain_slot: internal coding error");
        }
        if (current_domain.type() != TILEDB_NDRECTANGLE) {
            throw TileDBSOMAError(
                "_core_current_domain_slot: found non-rectangle type");
        }
        NDRectangle ndrect = current_domain.ndrectangle();

        // Convert from two-element array (core API) to pair (tiledbsoma API)
        std::array<std::string, 2> arr = ndrect.range<std::string>(name);

        // Here is an intersection of a few oddities:
        //
        // * Core domain for string dims must be a nullptr pair; it cannot be
        //   anything else.
        // * TileDB-Py shows this by using an empty-string pair, which we
        //   imitate.
        // * Core current domain for string dims must _not_ be a nullptr pair.
        // * In TileDB-SOMA, unless the user specifies otherwise, we use "" for
        //   min and "\xff" for max.
        // * However, "\xff" causes display problems in Python. It's also
        //   flat-out confusing to show to users.
        //
        // To work with all these factors, if the current domain is the default
        // "" to "\xff", return an empty-string pair just as we do for domain.
        if (arr[0] == "" && arr[1] == "\xff") {
            return std::pair<std::string, std::string>("", "");
        } else {
            return std::pair<std::string, std::string>(arr[0], arr[1]);
        }
    }

    /**
     * Returns the core domain at the given dimension.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     *   - soma maxdomain is core domain
     *   - core current domain is not accessed at the soma level
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    template <typename T>
    std::pair<T, T> _core_domain_slot(const std::string& name) const {
        if (std::is_same_v<T, std::string>) {
            throw std::runtime_error(
                "SOMAArray::_core_domain_slot: template-specialization "
                "failure.");
        }
        return arr_->schema().domain().dimension(name).domain<T>();
    }

    std::pair<std::string, std::string> _core_domain_slot_string(
        const std::string&) const {
        // Core domain for string dims is always a nullptr pair at the C++
        // level.  We follow the convention started by TileDB-Py which is to
        // report these as an empty-string pair.
        return std::pair<std::string, std::string>("", "");
    }

    /**
     * Returns the SOMA domain at the given dimension.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     */
    template <typename T>
    std::pair<T, T> soma_domain_slot(const std::string& name) const {
        if (has_current_domain()) {
            return _core_current_domain_slot<T>(name);
        } else {
            return _core_domain_slot<T>(name);
        }
    }

    /**
     * Returns the SOMA maxdomain at the given dimension.
     *
     * o For arrays with core current-domain support:
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma maxdomain is core domain
     */
    template <typename T>
    std::pair<T, T> soma_maxdomain_slot(const std::string& name) const {
        return _core_domain_slot<T>(name);
    }

    /**
     * Returns the SOMA domain in its entirety, as an Arrow table for return to
     * Python/R.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     *   - soma maxdomain is core domain
     *   - core current domain is not accessed at the soma level
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    ArrowTable get_soma_domain() {
        if (has_current_domain()) {
            return _get_core_current_domain();
        } else {
            return _get_core_domain();
        }
    }

    /**
     * Returns the SOMA maxdomain in its entirety, as an Arrow table for return
     * to Python/R.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     *   - soma maxdomain is core domain
     *   - core current domain is not accessed at the soma level
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    ArrowTable get_soma_maxdomain() {
        return _get_core_domain();
    }

    /**
     * Returns the core non-empty domain in its entirety, as an Arrow
     * table for return to Python/R.
     */
    ArrowTable get_non_empty_domain() {
        return _get_core_domainish(Domainish::kind_non_empty_domain);
    }

    /**
     * Code-dedupe helper for core domain, core current domain, and core
     * non-empty domain.
     */
    ArrowTable _get_core_domainish(enum Domainish which_kind);

    /**
     * This enables some code deduplication between core domain, core current
     * domain, and core non-empty domain.
     */
    template <typename T>
    std::pair<T, T> _core_domainish_slot(
        const std::string& name, enum Domainish which_kind) const {
        if (std::is_same_v<T, std::string>) {
            throw std::runtime_error(
                "SOMAArray::_core_domainish_slot: template-specialization "
                "failure.");
        }
        switch (which_kind) {
            case Domainish::kind_core_domain:
                return _core_domain_slot<T>(name);
            case Domainish::kind_core_current_domain:
                return _core_current_domain_slot<T>(name);
            case Domainish::kind_non_empty_domain:
                return non_empty_domain_slot<T>(name);
            default:
                throw std::runtime_error(
                    "internal coding error in SOMAArray::_core_domainish_slot: "
                    "unknown kind");
        }
    }

    std::pair<std::string, std::string> _core_domainish_slot_string(
        const std::string& name, enum Domainish which_kind) const {
        switch (which_kind) {
            case Domainish::kind_core_domain:
                return _core_domain_slot_string(name);
            case Domainish::kind_core_current_domain:
                return _core_current_domain_slot_string(name);
            case Domainish::kind_non_empty_domain:
                return non_empty_domain_slot_var(name);
            default:
                throw std::runtime_error(
                    "internal coding error in "
                    "SOMAArray::_core_domainish_slot_string: unknown kind");
        }
    }

    /**
     * @brief Get the total number of unique cells in the array.
     *
     * @return uint64_t Total number of unique cells
     */
    uint64_t nnz();

    /**
     * @brief Get the current capacity of each dimension.
     *
     * This applies to arrays all of whose dims are of type int64_t: this
     * includes SOMASparseNDArray and SOMADenseNDArray, and default-indexed
     * SOMADataFrame.
     *
     * At the TileDB-SOMA level we call this "shape". At the TileDB Core
     * storage level this maps to "current domain".
     *
     * Further, we map this single n to the pair (0, n-1) since core permits a
     * doubly inclusive pair (lo, hi) on each dimension slot.
     *
     * @return A vector with length equal to the number of dimensions; each
     * value in the vector is the capacity of each dimension.
     */
    std::vector<int64_t> shape();

    /**
     * @brief Get the maximum resizable capacity of each dimension.
     *
     * This applies to arrays all of whose dims are of type int64_t: this
     * includes SOMASparseNDArray and SOMADenseNDArray, and default-indexed
     * SOMADataFrame.
     *
     * At the TileDB-SOMA level we call this "maxshape". At the TileDB Core
     * storage level this maps to "domain".
     *
     * Further, we map this single n to the pair (0, n-1) since core permits a
     * doubly inclusive pair (lo, hi) on each dimension slot.
     *
     * @return A vector with length equal to the number of dimensions; each
     * value in the vector is the maximum capacity of each dimension.
     */
    std::vector<int64_t> maxshape();

    /**
     * This wires up to Python/R to tell a user if they can call resize() on an
     * array without error. For single arrays, they could just call resize() and
     * take their chances -- but for experiment-level resize (e.g. append mode)
     * it's crucial that we provide a can-we-do-them-all pass through all arrays
     * in the experiment before attempting any of them.
     *
     * On failure, returns false and an error string suitable for showing
     * to the user; on success, returns true and the empty string.
     *
     * Failure reasons: the requested shape's dimension-count doesn't match the
     * arrays; the array doesn't have a shape set (they must call
     * upgrade_shape), or the requested shape doesn't fit within the array's
     * existing core domain.
     */
    std::pair<bool, std::string> can_resize(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages) {
        return _can_set_shape_helper(
            newshape, true, function_name_for_messages);
    }

    /**
     * This wires up to Python/R to tell a user if they can call
     * upgrade_shape() on an array without error. For single dataframes,
     * they could just call upgrade_shape() and take their chances -- but for
     * experiment-level resize (e.g. append mode) it's crucial that we provide a
     * can-we-do-them-all pass through all arrays in the experiment before
     * attempting any of them.
     *
     * On failure, returns false and an error string suitable for showing
     * to the user; on success, returns true and the empty string.
     *
     * Failure reasons: the requested shape's dimension-count doesn't match the
     * arrays; the array already has a shape set (they must call resize), the
     * requested shape doesn't fit within the array's existing core domain, or
     * the requested shape is a downsize of the array's existing core current
     * domain.
     */
    std::pair<bool, std::string> can_upgrade_shape(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages) {
        return _can_set_shape_helper(
            newshape, false, function_name_for_messages);
    }

    /**
     * This is similar to can_upgrade_shape, but it's a can-we call
     * for maybe_resize_soma_joinid.
     */
    std::pair<bool, std::string> can_resize_soma_joinid_shape(
        int64_t newshape, std::string function_name_for_messages);

    /**
     * @brief Resize the shape (what core calls "current domain") up to the
     * maxshape (what core calls "domain").
     *
     * This applies to arrays all of whose dims are of type int64_t: this
     * includes SOMASparseNDArray and SOMADenseNDArray, and default-indexed
     * SOMADataFrame.
     *
     * @return Nothing. Raises an exception if the resize would be a downsize,
     * which is not supported.
     */
    void resize(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages);

    /**
     * @brief Given an old-style array without current domain, sets its
     * current domain. This is applicable only to arrays having all dims
     * of int64 type. Namely, all SparseNDArray/DenseNDArray, and
     * default-indexed DataFrame.
     */
    void upgrade_shape(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages);

    /**
     * @brief Increases the tiledbsoma shape up to at most the maxshape,
     * resizing the soma_joinid dimension if it is a dimension.
     *
     * While SOMA SparseNDArray and DenseNDArray, along with default-indexed
     * DataFrame, have int64_t dims, non-default-indexed DataFrame objects need
     * not: it is only required that they have a dim _or_ an attr called
     * soma_joinid. If soma_joinid is one of the dims, it will be resized while
     * the others will be preserved. If soma_joinid is not one of the dims,
     * nothing will be changed, as nothing _needs_ to be changed.
     *
     * @return Throws if the requested shape exceeds the array's create-time
     * maxshape. Throws if the array does not have current-domain support.
     */
    void resize_soma_joinid_shape(
        int64_t newshape, std::string function_name_for_messages);

   protected:
    // These two are for use nominally by SOMADataFrame. This could be moved in
    // its entirety to SOMADataFrame, but it would entail moving several
    // SOMAArray attributes from private to protected, which has knock-on
    // effects on the order of constructor initializers, etc.: in total it's
    // simplest to place this here and have SOMADataFrame invoke it.
    //
    // They return the shape and maxshape slots for the soma_joinid dim, if
    // the array has one. These are important test-points and dev-internal
    // access-points, in particular, for the tiledbsoma-io experiment-level
    // resizer.
    std::optional<int64_t> _maybe_soma_joinid_shape();
    std::optional<int64_t> _maybe_soma_joinid_maxshape();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    uint64_t _get_max_capacity(tiledb_datatype_t index_type);

    /**
     * Convenience function for creating an ArraySchemaEvolution object
     * referencing this array's context pointer, along with its open-at
     * timestamp (if any).
     */
    ArraySchemaEvolution _make_se();

    /**
     * The caller must check the return value for .is_empty() to see if this is
     * a new-style array with current-domain support (.is_empty() is false) , or
     * an old-style array without current-domain support (.is_empty() is true).
     * We could implement this as a std::optional<CurrentDomain> return value
     * here, but, that would be a redundant indicator.
     */
    CurrentDomain _get_current_domain() const {
        return tiledb::ArraySchemaExperimental::current_domain(
            *ctx_->tiledb_ctx(), arr_->schema());
    }

    /**
     * Returns the core current domain in its entirety, as an Arrow
     * table for return to Python/R.
     */
    ArrowTable _get_core_current_domain() {
        return _get_core_domainish(Domainish::kind_core_current_domain);
    }

    /**
     * Returns the core domain in its entirety, as an Arrow
     * table for return to Python/R.
     */
    ArrowTable _get_core_domain() {
        return _get_core_domainish(Domainish::kind_core_domain);
    }

    /**
     * This is a code-dedupe helper for can_resize and can_upgrade_domain.
     */
    std::pair<bool, std::string> _can_set_shape_helper(
        const std::vector<int64_t>& newshape,
        bool is_resize,
        std::string function_name_for_messages);

    /**
     * This is a second-level code-dedupe helper for _can_set_shape_helper.
     */
    std::pair<bool, std::string> _can_set_shape_domainish_subhelper(
        const std::vector<int64_t>& newshape,
        bool check_current_domain,
        std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for resize and upgrade_shape.
     */
    void _set_current_domain_from_shape(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages);

    /**
     * While SparseNDArray, DenseNDArray, and default-indexed DataFrame
     * have int64 dims, variant-indexed DataFrames do not. This helper
     * lets us pre-check any attempts to treat dims as if they were int64.
     */
    bool _dims_are_int64();

    /**
     * Same, but throws.
     */
    void _check_dims_are_int64();

    /**
     * With old shape: core domain used to map to tiledbsoma shape; core current
     * domain did not exist.
     *
     * With new shape: core domain maps to tiledbsoma maxshape;
     * core current_domain maps to tiledbsoma shape.
     *
     * Here we distinguish between user-side API, and core-side implementation.
     */
    std::vector<int64_t> _tiledb_domain();
    std::vector<int64_t> _tiledb_current_domain();
    std::optional<int64_t> _maybe_soma_joinid_tiledb_current_domain();
    std::optional<int64_t> _maybe_soma_joinid_tiledb_domain();

    bool _cast_column(
        ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se);

    void _promote_indexes_to_values(ArrowSchema* schema, ArrowArray* array);

    template <typename T>
    void _cast_dictionary_values(ArrowSchema* schema, ArrowArray* array);

    std::vector<int64_t> _get_index_vector(
        ArrowSchema* schema, ArrowArray* array) {
        auto index_type = ArrowAdapter::to_tiledb_format(schema->format);

        switch (index_type) {
            case TILEDB_INT8: {
                int8_t* idxbuf = (int8_t*)array->buffers[1];
                return std::vector<int64_t>(idxbuf, idxbuf + array->length);
            }
            case TILEDB_UINT8: {
                uint8_t* idxbuf = (uint8_t*)array->buffers[1];
                return std::vector<int64_t>(idxbuf, idxbuf + array->length);
            }
            case TILEDB_INT16: {
                int16_t* idxbuf = (int16_t*)array->buffers[1];
                return std::vector<int64_t>(idxbuf, idxbuf + array->length);
            }
            case TILEDB_UINT16: {
                uint16_t* idxbuf = (uint16_t*)array->buffers[1];
                return std::vector<int64_t>(idxbuf, idxbuf + array->length);
            }
            case TILEDB_INT32: {
                int32_t* idxbuf = (int32_t*)array->buffers[1];
                return std::vector<int64_t>(idxbuf, idxbuf + array->length);
            }
            case TILEDB_UINT32: {
                uint32_t* idxbuf = (uint32_t*)array->buffers[1];
                return std::vector<int64_t>(idxbuf, idxbuf + array->length);
            }
            case TILEDB_INT64: {
                int64_t* idxbuf = (int64_t*)array->buffers[1];
                return std::vector<int64_t>(idxbuf, idxbuf + array->length);
            }
            case TILEDB_UINT64: {
                uint64_t* idxbuf = (uint64_t*)array->buffers[1];
                return std::vector<int64_t>(idxbuf, idxbuf + array->length);
            }
            default:
                throw TileDBSOMAError(
                    "Saw invalid index type when trying to promote indexes to "
                    "values");
        }
    }

    template <typename UserType>
    bool _cast_column_aux(
        ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se);

    template <typename UserType, typename DiskType>
    bool _set_column(
        ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se) {
        // Here we cast the passed-in column to be what is the type in the
        // schema on disk. For columns with dictionaries, we will need
        // additional processing steps

        UserType* buf;
        if (array->n_buffers == 3) {
            buf = (UserType*)array->buffers[2] + array->offset;
        } else {
            buf = (UserType*)array->buffers[1] + array->offset;
        }

        bool has_attr = tiledb_schema()->has_attribute(schema->name);
        if (has_attr && attr_has_enum(schema->name)) {
            // For columns with dictionaries, we need to set the data buffers to
            // the dictionary's indexes. If there were any new enumeration
            // values added, we need to extend and and evolve the TileDB
            // ArraySchema

            // Return whether we extended the enumeration for this attribute
            return _extend_enumeration(
                schema->dictionary,  // value schema
                array->dictionary,   // value array
                schema,              // index schema
                array,               // index array
                se);
        } else {
            // In the general case, we can just cast the values and set the
            // write buffers
            std::vector<UserType> orig_vals(buf, buf + array->length);
            std::vector<DiskType> casted_values(
                orig_vals.begin(), orig_vals.end());

            mq_->setup_write_column(
                schema->name,
                casted_values.size(),
                (const void*)casted_values.data(),
                (uint64_t*)nullptr,
                (uint8_t*)array->buffers[0]);

            // Return false because we do not extend the enumeration
            return false;
        }
    }

    template <typename ValueType>
    bool _extend_and_evolve_schema(
        ArrowSchema* value_schema,
        ArrowArray* value_array,
        ArrowSchema* index_schema,
        ArrowArray* index_array,
        ArraySchemaEvolution se);

    template <typename ValueType>
    void _remap_indexes(
        std::string name,
        Enumeration extended_enmr,
        std::vector<ValueType> enums_in_write,
        ArrowSchema* index_schema,
        ArrowArray* index_array) {
        // If the passed-in enumerations are only a subset of the new extended
        // enumerations, then we will need to remap the indexes. Here identify
        // the dictionary values' type

        auto user_index_type = ArrowAdapter::to_tiledb_format(
            index_schema->format);
        switch (user_index_type) {
            case TILEDB_INT8:
                return _remap_indexes_aux<ValueType, int8_t>(
                    name, extended_enmr, enums_in_write, index_array);
            case TILEDB_UINT8:
                return _remap_indexes_aux<ValueType, uint8_t>(
                    name, extended_enmr, enums_in_write, index_array);
            case TILEDB_INT16:
                return _remap_indexes_aux<ValueType, int16_t>(
                    name, extended_enmr, enums_in_write, index_array);
            case TILEDB_UINT16:
                return _remap_indexes_aux<ValueType, uint16_t>(
                    name, extended_enmr, enums_in_write, index_array);
            case TILEDB_INT32:
                return _remap_indexes_aux<ValueType, int32_t>(
                    name, extended_enmr, enums_in_write, index_array);
            case TILEDB_UINT32:
                return _remap_indexes_aux<ValueType, uint32_t>(
                    name, extended_enmr, enums_in_write, index_array);
            case TILEDB_INT64:
                return _remap_indexes_aux<ValueType, int64_t>(
                    name, extended_enmr, enums_in_write, index_array);
            case TILEDB_UINT64:
                return _remap_indexes_aux<ValueType, uint64_t>(
                    name, extended_enmr, enums_in_write, index_array);
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
        // Get the user passed-in dictionary indexes
        IndexType* idxbuf;
        if (index_array->n_buffers == 3) {
            idxbuf = (IndexType*)index_array->buffers[2];
        } else {
            idxbuf = (IndexType*)index_array->buffers[1];
        }
        std::vector<IndexType> original_indexes(
            idxbuf, idxbuf + index_array->length);

        // Shift the dictionary indexes to match the on-disk extended
        // enumerations
        std::vector<IndexType> shifted_indexes;
        auto enmr_vec = extended_enmr.as_vector<ValueType>();
        for (auto i : original_indexes) {
            // For nullable columns, when the value is NULL, the associated
            // index may be a negative integer, so do not index into
            // enums_in_write or it will segfault
            if (0 > i) {
                shifted_indexes.push_back(i);
            } else {
                auto it = std::find(
                    enmr_vec.begin(), enmr_vec.end(), enums_in_write[i]);
                shifted_indexes.push_back(it - enmr_vec.begin());
            }
        }

        // Cast the user passed-in index type to be what is on-disk before we
        // set the write buffers. Here we identify the on-disk type
        auto attr = tiledb_schema()->attribute(column_name);
        switch (attr.type()) {
            case TILEDB_INT8:
                return _cast_shifted_indexes<IndexType, int8_t>(
                    column_name, shifted_indexes, index_array);
            case TILEDB_UINT8:
                return _cast_shifted_indexes<IndexType, uint8_t>(
                    column_name, shifted_indexes, index_array);
            case TILEDB_INT16:
                return _cast_shifted_indexes<IndexType, int16_t>(
                    column_name, shifted_indexes, index_array);
            case TILEDB_UINT16:
                return _cast_shifted_indexes<IndexType, uint16_t>(
                    column_name, shifted_indexes, index_array);
            case TILEDB_INT32:
                return _cast_shifted_indexes<IndexType, int32_t>(
                    column_name, shifted_indexes, index_array);
            case TILEDB_UINT32:
                return _cast_shifted_indexes<IndexType, uint32_t>(
                    column_name, shifted_indexes, index_array);
            case TILEDB_INT64:
                return _cast_shifted_indexes<IndexType, int64_t>(
                    column_name, shifted_indexes, index_array);
            case TILEDB_UINT64:
                return _cast_shifted_indexes<IndexType, uint64_t>(
                    column_name, shifted_indexes, index_array);
            default:
                throw TileDBSOMAError(
                    "Saw invalid enumeration index type when trying to extend"
                    "enumeration");
        }
    }

    template <typename UserIndexType, typename DiskIndexType>
    void _cast_shifted_indexes(
        std::string column_name,
        std::vector<UserIndexType> shifted_indexes,
        ArrowArray* index_array) {
        // Cast the user passed-in index type to be what is on-disk and
        // set the write buffers
        std::vector<DiskIndexType> casted_indexes(
            shifted_indexes.begin(), shifted_indexes.end());

        mq_->setup_write_column(
            column_name,
            casted_indexes.size(),
            (const void*)casted_indexes.data(),
            (uint64_t*)nullptr,
            (uint8_t*)index_array->buffers[0]);
    }

    bool _extend_enumeration(
        ArrowSchema* value_schema,
        ArrowArray* value_array,
        ArrowSchema* index_schema,
        ArrowArray* index_array,
        ArraySchemaEvolution se);

    void fill_metadata_cache();

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

    // Array associated with metadata_. Metadata values need to be
    // accessible in write mode as well. We need to keep this read-mode
    // array alive in order for the metadata value pointers in the cache to
    // be accessible
    std::shared_ptr<Array> meta_cache_arr_;

    // True if this is the first call to read_next()
    bool first_read_next_ = true;

    // True if the query was submitted
    bool submitted_ = false;

    // Unoptimized method for computing nnz() (issue `count_cells` query)
    uint64_t _nnz_slow();
};

// These are all specializations to string/bool of various methods
// which require special handling for that type.
//
// Declaring them down here is a bit weird -- they're easy to miss
// on a read-through. However, we're in a bit of a bind regarding
// various compilers: if we do these specializations within the
// `class SOMAArray { ... }`, then one compiler errors if we do
// include `template <>`, while another errors if we don't.
// Doing it down here, no compiler complains.

template <>
void SOMAArray::_cast_dictionary_values<std::string>(
    ArrowSchema* schema, ArrowArray* array);

template <>
void SOMAArray::_cast_dictionary_values<bool>(
    ArrowSchema* schema, ArrowArray* array);

template <>
bool SOMAArray::_cast_column_aux<std::string>(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se);

template <>
bool SOMAArray::_cast_column_aux<bool>(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se);

template <>
bool SOMAArray::_extend_and_evolve_schema<std::string>(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    ArraySchemaEvolution se);

}  // namespace tiledbsoma

#endif  // SOMA_ARRAY
