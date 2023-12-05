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

#include <tiledb/tiledb>

#include "../utils/common.h"
#include "array_buffers.h"
#include "column_buffer.h"

namespace tiledbsoma {

using namespace tiledb;

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

    ManagedQuery() = delete;
    ManagedQuery(const ManagedQuery&) = delete;
    ManagedQuery(ManagedQuery&&) = default;
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
     * @param names Vector of column names
     * @param if_not_empty Prevent changing an "empty" selection of all columns
     */
    void select_columns(
        const std::vector<std::string>& names, bool if_not_empty = false);

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
        subarray_range_set_ = true;
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
        subarray_range_set_ = true;
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
        subarray_range_set_ = true;
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
        subarray_range_set_ = true;
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
     * @param layout A tiledb_layout_t constant
     */
    void set_layout(tiledb_layout_t layout) {
        query_->set_layout(layout);
    }

    /**
     * @brief Set column data for write query.
     *
     * @param column_name Column name
     * @param buff Buffer array pointer with elements of the column type.
     * @param nelements Number of array elements in buffer
     */
    void set_column_data(
        std::string column_name, std::shared_ptr<ColumnBuffer> column_buffer) {
        if (array_->schema().array_type() == TILEDB_SPARSE ||
            schema_->has_attribute(column_name)) {
            auto data = column_buffer->data<std::byte>();
            query_->set_data_buffer(
                column_name, (void*)data.data(), data.size_bytes());
        } else {
            switch (column_buffer->type()) {
                case TILEDB_STRING_ASCII:
                case TILEDB_STRING_UTF8:
                case TILEDB_CHAR:
                case TILEDB_BLOB:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<std::string>()[0],
                        column_buffer->data<std::string>()[1]);
                    break;
                case TILEDB_FLOAT32:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<float>()[0],
                        column_buffer->data<float>()[1]);
                    break;
                case TILEDB_FLOAT64:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<double>()[0],
                        column_buffer->data<double>()[1]);
                    break;
                case TILEDB_UINT8:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<uint8_t>()[0],
                        column_buffer->data<uint8_t>()[1]);
                    break;
                case TILEDB_INT8:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<int8_t>()[0],
                        column_buffer->data<int8_t>()[1]);
                    break;
                case TILEDB_UINT16:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<uint16_t>()[0],
                        column_buffer->data<uint16_t>()[1]);
                    break;
                case TILEDB_INT16:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<int16_t>()[0],
                        column_buffer->data<int16_t>()[1]);
                    break;
                case TILEDB_UINT32:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<uint32_t>()[0],
                        column_buffer->data<uint32_t>()[1]);
                    break;
                case TILEDB_INT32:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<int32_t>()[0],
                        column_buffer->data<int32_t>()[1]);
                    break;
                case TILEDB_UINT64:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<uint64_t>()[0],
                        column_buffer->data<uint64_t>()[1]);
                    break;
                case TILEDB_INT64:
                case TILEDB_TIME_SEC:
                case TILEDB_TIME_MS:
                case TILEDB_TIME_US:
                case TILEDB_TIME_NS:
                case TILEDB_DATETIME_SEC:
                case TILEDB_DATETIME_MS:
                case TILEDB_DATETIME_US:
                case TILEDB_DATETIME_NS:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<int64_t>()[0],
                        column_buffer->data<int64_t>()[1]);
                    break;
                default:
                    break;
            }
            query_->set_subarray(*subarray_);
        }
    }

    /**
     * @brief Configure query and allocate result buffers for reads.
     *
     */
    void setup_read();

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
     * @brief Submit the write query.
     *
     */
    void submit_write();

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
        bool has_empty = false;
        for (auto subdim : subarray_range_empty_) {
            if (subdim.second == true) {
                has_empty = true;
                break;
            }
        }
        return subarray_range_set_ && has_empty;
    }

    /**
     * @brief Return the query type.
     *
     * @return TILEDB_READ or TILEDB_WRITE
     */
    tiledb_query_type_t query_type() const {
        return query_->query_type();
    }

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    /**
     * @brief Check if column name is contained in the query results.
     *
     * @param name Column name
     */
    void check_column_name(const std::string& name);

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
    bool subarray_range_set_ = false;

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
    std::future<void> query_future_;
};
};  // namespace tiledbsoma

#endif
