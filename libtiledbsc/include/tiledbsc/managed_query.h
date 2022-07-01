#ifndef MANAGED_QUERY_H
#define MANAGED_QUERY_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'
#include <unordered_set>

#include <tiledb/tiledb>

#include "tiledbsc/column_buffer.h"
#include "tiledbsc/common.h"

namespace tiledbsc {

using namespace tiledb;

constexpr size_t TILEDBSC_DEFAULT_ALLOC = 524288;

class ManagedQuery {
   public:
    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new Managed Query object
     *
     * @param array TileDB array
     * @param initial_cells Initial number of cells to allocate
     */
    ManagedQuery(
        std::shared_ptr<Array> array,
        size_t initial_cells = TILEDBSC_DEFAULT_ALLOC);

    /**
     * @brief Select columns names to query (dim and attr).
     *
     * @param names Vector of column names
     */
    void select_columns(std::vector<std::string> names);

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
            query_->add_range(dim, start, stop);
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
    void select_points(const std::string& dim, std::vector<T> points) {
        for (auto& point : points) {
            query_->add_range(dim, point, point);
        }
    }

    /**
     * @brief Execute the query and return the number of cells read.
     * To handle possible incomplete queries, `execute()` must be called until
     * it returns 0, which indicated the query is complete.
     *
     * For example:
     *   while (auto num_cells = mq.execute()) {
     *     // process the results
     *   }
     *
     * @return size_t Number of cells read. Returns 0 when the read is complete.
     */
    size_t execute();

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

    // Initial number of cells to allocate for each ColumnBuffer.
    size_t initial_cells_;

    // TileDB query being managed.
    std::unique_ptr<Query> query_;

    // Set of column names to read (dim and attr). If empty, query all columns.
    std::unordered_set<std::string> columns_;

    // Map of column name to ColumnBuffer.
    std::unordered_map<std::string, std::shared_ptr<ColumnBuffer>> buffers_;
};

};  // namespace tiledbsc

#endif
