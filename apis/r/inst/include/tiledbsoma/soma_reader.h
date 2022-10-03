/**
 * @file   soma_reader.h
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
 *   This declares the SOMAReader
 */

#ifndef SOMA_READER
#define SOMA_READER

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <future>

#include <tiledb/tiledb>

#include "thread_pool/thread_pool.h"
#include "tiledbsoma/managed_query.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAReader {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open an array at the specified URI and return SOMAReader object.
     *
     * @param uri URI of the array
     * @param name Name of the array
     * @param batch_size Read batch size
     * @param result_order Read result order
     * @param platform_config Config parameter dictionary
     * @return std::unique_ptr<SOMAReader> SOMAReader
     */
    static std::unique_ptr<SOMAReader> open(
        std::string_view uri,
        std::string_view name = "unnamed",
        std::map<std::string, std::string> platform_config = {},
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        std::string_view result_order = "auto");

    //===================================================================
    //= public non-static
    //===================================================================
    /**
     * @brief Construct a new SOMAReader object
     *
     * @param uri URI of the array
     * @param name name of the array
     * @param ctx TileDB context
     */
    SOMAReader(
        std::string_view uri,
        std::string_view name,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        std::string_view result_order);

    SOMAReader() = delete;
    SOMAReader(const SOMAReader&) = delete;
    SOMAReader(SOMAReader&&) = default;
    ~SOMAReader() = default;

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
        const std::span<T> points,
        int partition_index,
        int partition_count) {
        // Validate partition inputs
        if (partition_index >= partition_count) {
            throw TileDBSOMAError(fmt::format(
                "[SOMAReader] partition_index ({}) must be < partition_count "
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
                "[SOMAReader] set_dim_points partitioning: dim={} index={} "
                "count={} "
                "range=[{}, {}] of {} points",
                dim,
                partition_index,
                partition_count,
                start,
                start + partition_size - 1,
                points.size()));

            mq_->select_points(
                dim, std::span<T>{&points[start], partition_size});
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
     *   auto reader = SOMAReader::open(uri);
     *   reader->submit();
     *   while (auto batch = x_data->read_next()) {
     *       ...process batch ...
     *   }
     *
     * @return std::optional<std::shared_ptr<ArrayBuffers>>
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
        return mq_->results_complete();
    }

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMAReader URI
    std::string uri_;

    // Batch size
    std::string batch_size_;

    // Managed query for the array
    std::unique_ptr<ManagedQuery> mq_;

    // True if this is the first call to read_next()
    bool first_read_next_ = true;
};

}  // namespace tiledbsoma

#endif
