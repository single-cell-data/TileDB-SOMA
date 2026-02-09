/**
 * @file   managed_query.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This declares the managed query API.
 */

#ifndef COMMON_MANAGED_QUERY_H
#define COMMON_MANAGED_QUERY_H

#include <future>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "../concepts.h"

#pragma region Forward declarations

struct ArrowArray;
struct ArrowSchema;

namespace tiledb {
class Array;
class ArraySchemaEvolution;
class Attribute;
class Context;
class CurrentDomain;
class Enumeration;
class Query;
class QueryCondition;
class Subarray;
}  // namespace tiledb

namespace tiledbsoma::common {
class ArrayBuffers;
class ColumnBuffer;
class Status;
class WriteColumnBuffer;
}  // namespace tiledbsoma::common

#pragma endregion

namespace tiledbsoma::common {

enum class StatusCode { OK, ERROR };

class Status {
   public:
    Status();
    Status(StatusCode code, std::string_view origin, std::string_view message);
    ~Status() = default;

    StatusCode code() const {
        return code_;
    }

    std::string origin() const {
        return origin_;
    }

    std::string message() const {
        return message_;
    }

   private:
    StatusCode code_;
    std::string origin_;
    std::string message_;
};

enum class ResultOrder { AUTOMATIC = 0, ROWMAJOR, COLMAJOR, UNORDERED, GLOBAL };

class ManagedQuery {
   public:
    /**
     * @brief Construct a new ManagedQuery object
     *
     * @param array TileDB array
     * @param name Name of the array
     */
    ManagedQuery(
        std::shared_ptr<tiledb::Array> array, std::shared_ptr<tiledb::Context> ctx, std::string_view name = "unnamed");

    /** No default constructor. */
    ManagedQuery() = delete;

    /** No copy constructor: use move constructor instead. */
    ManagedQuery(const ManagedQuery&) = delete;

    /**
     * Default move constructor.
     *
     * Note: Each member should have a non-deleted move constructor. If there is a need to
     * violate this an explicitly defined move constructor should be added.
     */
    ManagedQuery(ManagedQuery&& other) = default;

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
     * @brief Set query result order (layout).
     *
     * @param layout A ResultOrder constant
     */
    void set_layout(ResultOrder layout);

    /**
     * @brief Returns the result order set by the query.
     *
     * @return ResultOrder
     */
    ResultOrder result_order() const;

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
    void select_columns(std::span<const std::string> names, bool if_not_empty = false, bool replace = false);

    /**
     * @brief Reset column selection to none, aka "all".
     */
    void reset_columns();

    std::vector<std::string> column_names() const;

#pragma region Member functions

#pragma region Read functions

    /**
     * @brief Select dimension ranges to query.
     *
     * @tparam T Dimension type
     * @param dim Dimension name
     * @param ranges Vector of dimension ranges
     */
    template <typename T>
    void select_ranges(const std::string& dim, std::span<const std::pair<T, T>> ranges) {
        subarray_range_set_[dim] = true;
        subarray_range_empty_[dim] = true;
        for (auto& [start, stop] : ranges) {
            ManagedQuery::select_range(subarray_.get(), dim, start, stop);
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
    void select_points(const std::string& dim, std::span<const T> points) {
        subarray_range_set_[dim] = true;
        subarray_range_empty_[dim] = true;
        for (auto& point : points) {
            ManagedQuery::select_range(subarray_.get(), dim, point, point);
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
        ManagedQuery::select_range(subarray_.get(), dim, point, point);
        subarray_range_set_[dim] = true;
        subarray_range_empty_[dim] = false;
    }

    std::optional<std::shared_ptr<ArrayBuffers>> read_next();

    /**
     * @brief Return true if the query result buffers hold all results from the
     * query. The return value is false if the query was incomplete.
     *
     * @return true The buffers hold all results from the query.
     */
    bool results_complete() const;

    uint64_t total_num_cells() const;

#pragma endregion

#pragma region Write functions

    /**
     * @brief Set column data for write query.
     *
     * @param name Column name
     * @param num_elems Number of array elements in buffer
     * @param data Pointer to the data buffer
     *  If the data type is Boolean, the data has already been casted to uint8
     * @param offsets Pointer to the offsets buffer
     * @param validity Vector of validity buffer casted to uint8
     */
    template <typename DataStorage, typename OffsetStorage>
        requires is_data_buffer<DataStorage> && is_offset_buffer<OffsetStorage>
    void setup_write_column(
        std::string_view name,
        uint64_t num_elems,
        DataStorage data,
        OffsetStorage offsets,
        std::unique_ptr<uint8_t[]> validity = nullptr,
        bool copy_buffers = false) {
        // Create the ArrayBuffers as necessary
        if (!buffers_) {
            buffers_ = std::make_shared<ArrayBuffers>();
        }

        ManagedQuery::set_write_column(
            buffers_.get(),
            *array_,
            *query_,
            subarray_.get(),
            name,
            num_elems,
            std::move(data),
            std::move(offsets),
            std::move(validity),
            copy_buffers);
    }

    /**
     * @brief Set the write buffers for an Arrow Table or Batch as represented
     * by an ArrowSchema and ArrowArray. Nulls values are not allowed for
     * dimensions and will error out. Any null values present in non-nullable
     * attributes will be casted to fill values for the given TileDB datatype.
     *
     * @param arrow_schema
     * @param arrow_array
     */
    void set_array_data(ArrowSchema* arrow_schema, ArrowArray* arrow_array);

    /**
     * @brief Submit the write query.
     *
     */
    void submit_write();

    /**
     * @brief Finalize the write query.
     *
     */
    void finalize();

    /**
     * @brief Submit and finalize the write query.
     *
     */
    void submit_and_finalize();

    /**
     * @brief Set a query condition.
     *
     * @param qc A TileDB QueryCondition
     */
    void set_condition(const tiledb::QueryCondition& qc);

#pragma endregion

    /**
     * @brief Extend an Enumeration with the values contained in the supplied Arrow array and return true if the enumeration is modified
     */
    std::optional<std::shared_ptr<tiledb::Enumeration>> extend_enumeration(
        ArrowSchema* value_schema,
        ArrowArray* value_array,
        ArrowSchema* index_schema,
        ArrowArray* index_array,
        std::string_view column_name,
        std::string_view enumeration_name,
        tiledb::ArraySchemaEvolution& se,
        bool deduplicate);

    std::shared_ptr<tiledb::Array> array() const;

    std::shared_ptr<ArrayBuffers> buffers() const;

    /**
     * @brief Return true if the only ranges selected were empty.
     *
     * @return true if the query contains only empty ranges.
     */
    bool is_empty_query() const;

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
    bool is_complete(bool query_status_only = false) const;

#pragma endregion

    static const size_t MAX_RETRIES = 3;

   private:
#pragma region Type-erased static helpers

    template <typename T, std::same_as<tiledb::Subarray> Subarray = tiledb::Subarray>
    static void select_range(Subarray* subarray, const std::string& dim, T start, T stop) {
        subarray->add_range(dim, start, stop);
    }

    template <
        typename DataStorage,
        typename OffsetStorage,
        std::same_as<ArrayBuffers> Buffers = ArrayBuffers,
        std::same_as<WriteColumnBuffer> Column = WriteColumnBuffer>
        requires is_data_buffer<DataStorage> && is_offset_buffer<OffsetStorage>
    static void set_write_column(
        Buffers* buffers,
        const tiledb::Array& array,
        tiledb::Query& query,
        tiledb::Subarray* subarray,
        std::string_view name,
        uint64_t num_elems,
        DataStorage data,
        OffsetStorage offsets,
        std::unique_ptr<uint8_t[]> validity = nullptr,
        bool copy_buffers = false) {
        buffers->emplace(
            std::string(name),
            Column::create(
                array, name, num_elems, std::move(data), std::move(offsets), std::move(validity), copy_buffers));

        buffers->at(std::string(name))->attach(query, subarray);
    }

#pragma endregion

#pragma region Member functions

#pragma region Read functions

    /**
     * @brief Configure query and allocate result buffers for reads.
     */
    void setup_read();

    void submit_read();

#pragma endregion

#pragma region Write functions

    void setup_write();

    void teardown_write();

    bool setup_write_arrow_column(ArrowSchema* schema, ArrowArray* array, tiledb::ArraySchemaEvolution& se);

    void remap_enumeration_indices(
        ArrowSchema* schema,
        ArrowArray* array,
        const tiledb::Attribute& attribute,
        const tiledb::Enumeration& enumeration);

    void promote_indexes_to_values(ArrowSchema* schema, ArrowArray* array);

    /**
     * @brief Get the mapping of attributes to Enumerations.
     *
     * @return std::map<std::string, Enumeration>
     */
    std::map<std::string, tiledb::Enumeration> attribute_to_enumeration_mapping() const;

#pragma endregion

    void fill_in_subarrays();

    void fill_in_subarrays_domain();

    void fill_in_subarrays_current_domain(const tiledb::CurrentDomain& current_domain);

    bool has_any_subarray_range_set() const;

    bool has_any_empty_range() const;

#pragma endregion

#pragma region Member variables

    // TileDB context object
    std::shared_ptr<tiledb::Context> ctx_;

    // TileDB array being queried.
    std::shared_ptr<tiledb::Array> array_;

    // Name displayed in log messages
    std::string name_;

    // Set of column names to read (dim and attr). If empty, query all columns.
    std::vector<std::string> columns_;

    // A collection of ColumnBuffers attached to the query
    std::shared_ptr<ArrayBuffers> buffers_;

    // TileDB query being managed.
    std::unique_ptr<tiledb::Query> query_;

    // TileDB subarray containing the ranges for slicing.
    std::unique_ptr<tiledb::Subarray> subarray_;

    // True if a range has been added to the subarray
    std::map<std::string, bool> subarray_range_set_ = {};

    // Map whether the dimension is empty (true) or not
    std::map<std::string, bool> subarray_range_empty_ = {};

    // Future for asyncronous query
    std::future<Status> query_future_;

    // Query layout
    ResultOrder layout_ = ResultOrder::AUTOMATIC;

    // Number of query submission that returned no results due to insufficient buffer sizes
    size_t retries = 0;

    // Results in the buffers are complete (the query was never incomplete)
    bool results_complete_ = true;

    // Total number of cells read by the query
    uint64_t total_num_cells_ = 0;

    // True if the query has been submitted
    bool query_submitted_ = false;

#pragma endregion
};
}  // namespace tiledbsoma::common

#endif