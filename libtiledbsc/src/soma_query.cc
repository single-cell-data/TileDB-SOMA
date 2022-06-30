#include "tiledbsc/soma_query.h"
#include "tiledbsc/common.h"
#include "tiledbsc/logger_public.h"

namespace tiledbsc {
using namespace tiledb;

SOMAQuery::SOMAQuery(SOMA* soma, size_t index_alloc, size_t x_alloc) {
    mq_obs_ = std::make_unique<ManagedQuery>(
        soma->open_array("obs"), index_alloc);
    mq_var_ = std::make_unique<ManagedQuery>(
        soma->open_array("var"), index_alloc);
    mq_x_ = std::make_unique<ManagedQuery>(soma->open_array("X/data"), x_alloc);
}

std::optional<std::unordered_map<std::string, std::shared_ptr<ColumnBuffer>>>
SOMAQuery::next_results() {
    // Query is complete, return empty results
    if (mq_x_->status() == Query::Status::COMPLETE) {
        return std::nullopt;
    }

    // Abort if an invalid column was selected for the obs or var query.
    if (mq_obs_->is_invalid() || mq_var_->is_invalid()) {
        LOG_DEBUG(
            fmt::format("[SOMAQuery] Abort due to invalid selected column."));
        return std::nullopt;
    }

    if (mq_x_->status() == Query::Status::UNINITIALIZED) {
        // Submit obs and var query tasks in parallel
        auto obs_task = std::async(std::launch::async, [&]() {
            return query_and_select(mq_obs_, "obs_id");
        });

        auto var_task = std::async(std::launch::async, [&]() {
            return query_and_select(mq_var_, "var_id");
        });

        // Block until obs and var tasks complete
        auto num_obs = obs_task.get();
        auto num_var = var_task.get();

        // Return empty results if obs or var query was empty
        if (!num_obs || !num_var) {
            return std::nullopt;
        }
    }

    // Submit X query
    auto num_cells = mq_x_->submit();
    LOG_DEBUG(fmt::format("*** X cells read = {}", num_cells));

    // TODO: Build and return ArrowTable
    return mq_x_->results();
}

size_t SOMAQuery::query_and_select(
    std::unique_ptr<ManagedQuery>& mq, const std::string& dim_name) {
    // Select the dimension column, if all columns are not selected already.
    mq->select_columns({dim_name}, true);

    // Submit the query.
    auto num_cells = mq->submit();
    LOG_DEBUG(fmt::format("*** {} cells read = {}", dim_name, num_cells));

    // Throw an error if the query is incomplete.
    if (!mq->is_complete()) {
        throw TileDBSCError(fmt::format(
            "[SOMAQuery] {} query is incomplete. Increase buffer size.",
            dim_name));
    }

    // Return early if the query result is empty.
    if (num_cells == 0) {
        return num_cells;
    }

    // If the query contains a subset of the dimensions, apply the results to
    // the X query.
    // TODO: [optimize] If sliced query returns all results, skip this if block
    if (mq->is_sliced()) {
        auto points = mq->strings(dim_name);
        if (points.size() != num_cells) {
            throw TileDBSCError(fmt::format(
                "[SOMAQuery] {} query failed sanity check: {} != {}",
                dim_name,
                points.size(),
                num_cells));
        }

        // Add dimension range points to the X query.
        std::lock_guard<std::mutex> lock(mtx_);
        mq_x_->select_points(dim_name, points);
    }

    return num_cells;
}

}  // namespace tiledbsc
