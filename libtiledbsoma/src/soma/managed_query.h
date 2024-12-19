/**
 * @file   managed_query.h
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
 *   This declares the managed query API.
 */

#ifndef MANAGED_QUERY_H
#define MANAGED_QUERY_H

#include <future>
#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'
#include <unordered_set>
#include "span/span.hpp"

#include <tiledb/tiledb>

#include "../utils/common.h"
#include "array_buffers.h"
#include "column_buffer.h"
#include "enums.h"
#include "logger_public.h"

namespace tiledbsoma {

using namespace tiledb;
class SOMAArray;

// Probably we should just use a std::tuple here
class StatusAndException {
   public:
    StatusAndException(bool succeeded, std::string message)
        : succeeded_(succeeded)
        , message_(message) {
    }

    bool succeeded() {
        return succeeded_;
    }
    std::string message() {
        return message_;
    }

   private:
    bool succeeded_;
    std::string message_;
};

class ManagedQuery {
   public:
    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new ManagedQuery object
     *
     * @param array TileDB array
     * @param name Name of the array
     */
    ManagedQuery(
        std::shared_ptr<Array> array,
        std::shared_ptr<Context> ctx,
        std::string_view name = "unnamed");

    ManagedQuery(
        std::unique_ptr<SOMAArray> array,
        std::shared_ptr<Context> ctx,
        std::string_view name = "unnamed");

    ManagedQuery() = delete;

    ManagedQuery(const ManagedQuery&) = delete;

    ManagedQuery(ManagedQuery&& other)
        : array_(other.array_)
        , ctx_(other.ctx_)
        , name_(other.name_)
        , schema_(other.schema_)
        , query_(std::make_unique<Query>(*other.ctx_, *other.array_))
        , subarray_(std::make_unique<Subarray>(*other.ctx_, *other.array_))
        , subarray_range_set_(other.subarray_range_set_)
        , subarray_range_empty_(other.subarray_range_empty_)
        , columns_(other.columns_)
        , results_complete_(other.results_complete_)
        , total_num_cells_(other.total_num_cells_)
        , buffers_(other.buffers_)
        , query_submitted_(other.query_submitted_) {
    }

    ~ManagedQuery() = default;

    /**
     * @brief Close the array after waiting for any asynchronous queries to
     * complete.
     *
     */
    void close();

    /**
     * @brief Reset the state of this ManagedQuery object to prepare for a new
     * query, while holding the array open.
     *
     */
    void reset();

    /**
     * @brief Select columns names to query (dim and attr). If the
     * `if_not_empty` parameter is `true`, the column will be selected iff the
     * list of selected columns is empty. This prevents a `select_columns` call
     * from changing an empty list (all columns) to a subset of columns.
     *
     * NB: you may only select a given column once. Selecting twice will
     * generate an error in read_next.
     *
     * @param names Vector of column names
     * @param if_not_empty Prevent changing an "empty" selection of all columns
     * @param replace Column names will replace any existing selected columns.
     */
    void select_columns(
        const std::vector<std::string>& names,
        bool if_not_empty = false,
        bool replace = false);

    /**
     * @brief Reset column selection to none, aka "all".
     */
    void reset_columns();

    /**
     * @brief Returns the column names set by the query.
     *
     * @return std::vector<std::string>
     */
    std::vector<std::string> column_names() {
        return columns_;
    }

    /**
     * @brief Select dimension ranges to query.
     *
     * @tparam T Dimension type
     * @param dim Dimension name
     * @param ranges Vector of dimension ranges
     */
    template <typename T>
    void select_ranges(
        const std::string& dim, const std::vector<std::pair<T, T>>& ranges) {
        subarray_range_set_[dim] = true;
        subarray_range_empty_[dim] = true;
        for (auto& [start, stop] : ranges) {
            subarray_->add_range(dim, start, stop);
            subarray_range_empty_[dim] = false;
        }
    }

    /**
     * @brief Select dimension points to query.
     *
     * @tparam T Dimension type
     * @param dim Dimension name
     * @param points Vector of dimension points
     */
    template <typename T>
    void select_points(const std::string& dim, const std::vector<T>& points) {
        subarray_range_set_[dim] = true;
        subarray_range_empty_[dim] = true;
        for (auto& point : points) {
            subarray_->add_range(dim, point, point);
            subarray_range_empty_[dim] = false;
        }
    }

    /**
     * @brief Select dimension points to query.
     *
     * @tparam T Dimension type
     * @param dim Dimension name
     * @param points Vector of dimension points
     */
    template <typename T>
    void select_points(const std::string& dim, const tcb::span<T> points) {
        subarray_range_set_[dim] = true;
        subarray_range_empty_[dim] = true;
        for (auto& point : points) {
            subarray_->add_range(dim, point, point);
            subarray_range_empty_[dim] = false;
        }
    }

    /**
     * @brief Select dimension point to query.
     *
     * @tparam T Dimension type
     * @param dim Dimension name
     * @param point Dimension points
     */
    template <typename T>
    void select_point(const std::string& dim, const T& point) {
        subarray_->add_range(dim, point, point);
        subarray_range_set_[dim] = true;
        subarray_range_empty_[dim] = false;
    }

    /**
     * @brief Set a query condition.
     *
     * @param qc A TileDB QueryCondition
     */
    void set_condition(const QueryCondition& qc) {
        query_->set_condition(qc);
    }

    /**
     * @brief Set query result order (layout).
     *
     * @param layout A ResultOrder constant
     */
    void set_layout(ResultOrder layout);

    /**
     * @brief Set column data for write query.
     *
     * @param name Column name
     * @param num_elems Number of array elements in buffer
     * @param data Pointer to the data buffer
     * @param offsets Pointer to the offsets buffer
     * @param validity Pointer to the validity buffer
     */
    template <typename T>
    void setup_write_column(
        std::string_view name,
        uint64_t num_elems,
        const void* data,
        T* offsets,
        uint8_t* validity) {
        // Ensure the offset type is either uint32_t* or uint64_t*
        static_assert(
            std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>,
            "offsets must be either uint32_t* or uint64_t*");

        // Create the ArrayBuffers as necessary
        if (buffers_ == nullptr) {
            buffers_ = std::make_shared<ArrayBuffers>();
        }

        auto column = ColumnBuffer::create(array_, name);
        column->set_data(num_elems, data, offsets, validity);
        buffers_->emplace(std::string(name), column);
        buffers_->at(std::string(name))->attach(*query_, *subarray_);
    }

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
        return query_->query_status() == Query::Status::COMPLETE ||
               (!query_status_only && is_empty_query());
    }

    /**
     * @brief Return true if the query result buffers hold all results from the
     * query. The return value is false if the query was incomplete.
     *
     * @return true The buffers hold all results from the query.
     */
    bool results_complete() {
        return is_complete() && results_complete_;
    }

    /**
     * @brief Returns the total number of cells read so far, including any
     * previous incomplete queries.
     *
     * @return size_t Total number of cells read
     */
    size_t total_num_cells() {
        return total_num_cells_;
    }

    /**
     * @brief Return a view of data in column `name`.
     *
     * @tparam T Data type
     * @param name Column name
     * @return tcb::span<T> Data view
     */
    template <typename T>
    tcb::span<T> data(const std::string& name) {
        check_column_name(name);
        return buffers_->at(name)->data<T>();
    }

    /**
     * @brief Return a view of validity values for column `name`.
     *
     * @param name Column name
     * @return tcb::span<uint8_t> Validity view
     */
    const tcb::span<uint8_t> validity(const std::string& name) {
        check_column_name(name);
        return buffers_->at(name)->validity();
    }

    /**
     * @brief Return a vector of strings from the column `name`.
     *
     * @param name Column name
     * @return std::vector<std::string> Strings
     */
    std::vector<std::string> strings(const std::string& name) {
        check_column_name(name);
        return buffers_->at(name)->strings();
    }

    /**
     * @brief Return a string_view of the string at the provided cell index from
     * column `name`.
     *
     * @param name Column name
     * @param index Cell index
     * @return std::string_view String view
     */
    std::string_view string_view(const std::string& name, uint64_t index) {
        check_column_name(name);
        return buffers_->at(name)->string_view(index);
    }

    /**
     * @brief Submit the write query.
     *
     */
    void submit_write(bool sort_coords = true);

    /**
     * @brief Get the schema of the array.
     *
     * @return std::shared_ptr<ArraySchema> Schema
     */
    std::shared_ptr<ArraySchema> schema() {
        return schema_;
    }

    /**
     * @brief Return true if the only ranges selected were empty.
     *
     * @return true if the query contains only empty ranges.
     */
    bool is_empty_query() {
        return _has_any_empty_range() && _has_any_subarray_range_set();
    }

    /**
     * @brief Return the query type.
     *
     * @return TILEDB_READ or TILEDB_WRITE
     */
    tiledb_query_type_t query_type() const {
        return query_->query_type();
    }

    /**
     * @brief Return the query status.
     *
     * @return tiledb::Query::Status INCOMPLETE, COMPLETE, INPROGRESS, FAILED,
     * UNINITIALIZED, or INITIALIZED
     */
    Query::Status query_status() const {
        return query_->query_status();
    }

    bool is_first_read() const {
        return !query_submitted_;
    }

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    /**
     * @brief Configure query and allocate result buffers for reads.
     *
     */
    void setup_read();

    /**
     * @brief Submit the query.
     *
     */
    void submit_read();

    /**
     * @brief Return results from the query.
     *
     * @return std::shared_ptr<ArrayBuffers>
     */
    std::shared_ptr<ArrayBuffers> results();

    /**
     * @brief Check if column name is contained in the query results.
     *
     * @param name Column name
     */
    void check_column_name(const std::string& name);

    // Helper for is_empty_query
    bool _has_any_empty_range() {
        for (auto subdim : subarray_range_empty_) {
            if (subdim.second == true) {
                return true;
            }
        }
        return false;
    }

    // Helper for is_empty_query
    bool _has_any_subarray_range_set() {
        for (auto subdim : subarray_range_set_) {
            if (subdim.second == true) {
                return true;
            }
        }
        return false;
    }

    /**
     * This handles a few internals.
     *
     * One is that a dense array must have _at least one_
     * dim's subarray set for a read query. Without that, reads fail immediately
     * with the unambiguous
     *
     *   DenseReader: Cannot initialize reader; Dense reads must have a subarray
     *   set
     *
     * The other is a combination of several things. Firstly, is current-domain
     * support which we have for sparse arrays as of core 2.26, and for dense as
     * of 2.27. Secondly, without current-domain support, we had small domains;
     * with it, we have huge core domains (2^63-ish) which are immutable, and
     * small current domains which are upward-mutable. (The soma domain and
     * maxdomain, respectively, are core current domain and domain.) Thirdly,
     * if a query doesn't have a subarray set on any
     * particular dim, core will use the core domain on that dim. That was fine
     * when core domains were small; not fine now that they are huge. In this
     * routine, if the array is dense, for each dim without a subarray set,
     * we set it to match the soma domain. This guarantees correct behavior.
     */
    void _fill_in_subarrays_if_dense(bool is_read);
    void _fill_in_subarrays_if_dense_with_new_shape(
        const CurrentDomain& current_domain, bool is_read);
    void _fill_in_subarrays_if_dense_without_new_shape(bool is_read);

    // TileDB array being queried.
    std::shared_ptr<Array> array_;

    // TileDB context object
    std::shared_ptr<Context> ctx_;

    // Name displayed in log messages
    std::string name_;

    // Array schema
    std::shared_ptr<ArraySchema> schema_;

    // TileDB query being managed.
    std::unique_ptr<Query> query_;

    // TileDB subarray containing the ranges for slicing.
    std::unique_ptr<Subarray> subarray_;

    // True if a range has been added to the subarray
    std::map<std::string, bool> subarray_range_set_ = {};

    // Map whether the dimension is empty (true) or not
    std::map<std::string, bool> subarray_range_empty_ = {};

    // Set of column names to read (dim and attr). If empty, query all columns.
    std::vector<std::string> columns_;

    // Results in the buffers are complete (the query was never incomplete)
    bool results_complete_ = true;

    // Total number of cells read by the query
    size_t total_num_cells_ = 0;

    // A collection of ColumnBuffers attached to the query
    std::shared_ptr<ArrayBuffers> buffers_;

    // True if the query has been submitted
    bool query_submitted_ = false;

    // Future for asyncronous query
    std::future<StatusAndException> query_future_;

    /**
     * Convenience function for creating an ArraySchemaEvolution object
     * referencing this array's context pointer, along with its open-at
     * timestamp (if any).
     */
    ArraySchemaEvolution _make_se();

    uint64_t _get_max_capacity(tiledb_datatype_t index_type);

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
        // The array->offset is non-zero when we are passed sliced
        // Arrow tables like arrow_table[:m] or arrow_table[m:].
        if (array->n_buffers == 3) {
            buf = (UserType*)array->buffers[2] + array->offset;
        } else {
            buf = (UserType*)array->buffers[1] + array->offset;
        }
        uint8_t* validity = (uint8_t*)array->buffers[0];
        if (validity != nullptr) {
            validity += array->offset;
        }

        bool has_attr = schema_->has_attribute(schema->name);
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

            setup_write_column(
                schema->name,
                casted_values.size(),
                (const void*)casted_values.data(),
                (uint64_t*)nullptr,
                validity);

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
        auto attr = schema_->attribute(column_name);
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

        setup_write_column(
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

    /**
     * @brief Get the mapping of attributes to Enumerations.
     *
     * @return std::map<std::string, Enumeration>
     */
    std::map<std::string, Enumeration> get_attr_to_enum_mapping() {
        std::map<std::string, Enumeration> result;
        for (uint32_t i = 0; i < schema_->attribute_num(); ++i) {
            auto attr = schema_->attribute(i);
            if (attr_has_enum(attr.name())) {
                auto enmr_label = *get_enum_label_on_attr(attr.name());
                auto enmr = ArrayExperimental::get_enumeration(
                    *ctx_, *array_, enmr_label);
                result.insert({attr.name(), enmr});
            }
        }
        return result;
    }

    /**
     * @brief Get the Enumeration name associated with the given Attr.
     *
     * @return std::optional<std::string> The enumeration name if one exists.
     */
    std::optional<std::string> get_enum_label_on_attr(std::string attr_name) {
        auto attr = schema_->attribute(attr_name);
        return AttributeExperimental::get_enumeration_name(*ctx_, attr);
    }

    /**
     * @brief Check if the given attribute has an associated enumeration.
     *
     * @return bool
     */
    bool attr_has_enum(std::string attr_name) {
        return get_enum_label_on_attr(attr_name).has_value();
    }
};

// These are all specializations to string/bool of various methods
// which require special handling for that type.
//
// Declaring them down here is a bit weird -- they're easy to miss
// on a read-through. However, we're in a bit of a bind regarding
// various compilers: if we do these specializations within the
// `class ManagedQuery { ... }`, then one compiler errors if we do
// include `template <>`, while another errors if we don't.
// Doing it down here, no compiler complains.

template <>
void ManagedQuery::_cast_dictionary_values<std::string>(
    ArrowSchema* schema, ArrowArray* array);

template <>
void ManagedQuery::_cast_dictionary_values<bool>(
    ArrowSchema* schema, ArrowArray* array);

template <>
bool ManagedQuery::_cast_column_aux<std::string>(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se);

template <>
bool ManagedQuery::_cast_column_aux<bool>(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se);

template <>
bool ManagedQuery::_extend_and_evolve_schema<std::string>(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    ArraySchemaEvolution se);

};  // namespace tiledbsoma

#endif