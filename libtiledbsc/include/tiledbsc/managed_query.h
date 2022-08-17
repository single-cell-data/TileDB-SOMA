#ifndef MANAGED_QUERY_H
#define MANAGED_QUERY_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'
#include <unordered_set>

#include <tiledb/tiledb>

#include "tiledbsc/column_buffer.h"
#include "tiledbsc/common.h"

namespace tiledbsc {

using namespace tiledb;

class ManagedQuery {
   public:
    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new Managed Query object
     *
     * @param array TileDB array
     * @param name Name of the array
     */
    ManagedQuery(std::shared_ptr<Array> array, std::string name = "array");

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
        std::vector<std::string> names, bool if_not_empty = false);

    /**
     * @brief Select dimension ranges to query.
     *
     * @tparam T Dimension type
     * @param dim Dimension name
     * @param ranges Vector of dimension ranges
     */
    template <typename T>
    void select_ranges(
        const std::string& dim, std::vector<std::pair<T, T>> ranges) {
        for (auto& [start, stop] : ranges) {
            subarray_->add_range(dim, start, stop);
        }
        sliced_ = true;
    }

    /**
     * @brief Select dimension points to query.
     *
     * @tparam T Dimension type
     * @param dim Dimension name
     * @param points Vector of dimension points
     */
    template <typename T>
    void select_points(const std::string& dim, std::vector<T> points) {
        for (auto& point : points) {
            subarray_->add_range(dim, point, point);
        }
        sliced_ = true;
    }

    /**
     * @brief Set a query condition.
     *
     * @param qc A TileDB QueryCondition
     */
    void set_condition(const QueryCondition& qc) {
        query_->set_condition(qc);
        sliced_ = true;
    }

    /**
     * @brief Submit the query and return the number of cells read.
     * To handle incomplete queries, `submit()` must be called until
     * `is_complete()` is true.
     *
     * For example:
     *   while (!mq.is_complete()) {
     *     auto num_cells = mq.submit();
     *     // process results
     *   }
     *
     * @return size_t Number of cells read. Returns 0 when the read is complete.
     */
    size_t submit();

    /**
     * @brief Return the query status.
     *
     * @return Query::Status Query status
     */
    Query::Status status() {
        return query_->query_status();
    }

    /**
     * @brief Check if the query is complete.
     *
     * @return true Query status is COMPLETE
     */
    bool is_complete() {
        return query_->query_status() == Query::Status::COMPLETE;
    }

    /**
     * @brief Return true if an invalid column has been selected.
     *
     * @return true An invalid column was selected
     */
    bool is_invalid() {
        return invalid_columns_selected_;
    }

    /**
     * @brief Return true if the query has dimension ranges selected or a query
     * condition applied.
     *
     * @return true The query will be sliced
     */
    bool is_sliced() {
        return sliced_;
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
     * @return std::span<T> Data view
     */
    template <typename T>
    std::span<T> data(const std::string& name) {
        check_column_name(name);
        return buffers_.at(name)->data<T>();
    }

    /**
     * @brief Return a vector of strings from the column `name`.
     *
     * @param name Column name
     * @return std::vector<std::string> Strings
     */
    std::vector<std::string> strings(const std::string& name) {
        check_column_name(name);
        return buffers_.at(name)->strings();
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
        return buffers_.at(name)->string_view(index);
    }

    /**
     * @brief Return results from the query.
     *
     * @return ArrayBuffers Results
     */
    ArrayBuffers results() {
        return buffers_;
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
    void check_column_name(const std::string& name) {
        if (!buffers_.contains(name)) {
            throw TileDBSCError(fmt::format(
                "[ManagedQuery] Column '{}' is not available in the query "
                "results.",
                name));
        }
    }

    // TileDB array being queried.
    std::shared_ptr<Array> array_;

    // Name displayed in log messages
    std::string name_;

    // Array schema
    ArraySchema schema_;

    // TileDB query being managed.
    std::unique_ptr<Query> query_;

    // TileDB subarray containing the ranges for slicing.
    std::unique_ptr<Subarray> subarray_;

    // Set of column names to read (dim and attr). If empty, query all columns.
    std::unordered_set<std::string> columns_;

    // Invalid columns have been selected.
    bool invalid_columns_selected_ = false;

    // Query is sliced with dimension range slices or query conditions.
    bool sliced_ = false;

    // Results in the buffers are complete (the query was never incomplete)
    bool results_complete_ = true;

    // Total number of cells read by the query
    size_t total_num_cells_ = 0;

    // Map of column name to ColumnBuffer.
    ArrayBuffers buffers_;
};

};  // namespace tiledbsc

#endif
