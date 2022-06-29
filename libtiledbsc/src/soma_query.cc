#include "tiledbsc/soma_query.h"
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
    if (mq_x_->status() == Query::Status::COMPLETE) {
        return std::nullopt;
    }

    if (mq_x_->status() == Query::Status::UNINITIALIZED) {
        auto num_obs = query_and_select(mq_obs_, "obs_id");
        auto num_var = query_and_select(mq_var_, "var_id");

        // Return empty results if obs or var query was empty
        if (!num_obs || !num_var) {
            return std::nullopt;
        }
    }

    auto num_cells = mq_x_->submit();
    LOG_DEBUG(fmt::format("*** X cells read = {}", num_cells));

    std::unordered_map<std::string, std::shared_ptr<ColumnBuffer>> results;

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

    auto points = mq->strings(dim_name);
    if (points.size() != num_cells) {
        throw TileDBSCError(fmt::format(
            "[SOMAQuery] {} query failed sanity check: {} != {}",
            dim_name,
            points.size(),
            num_cells));
    }

    // Add dimension range points to the X query.
    std::lock_guard<std::mutex> lck(mtx_);
    mq_x_->select_points(dim_name, points);

    return num_cells;
}

}  // namespace tiledbsc
