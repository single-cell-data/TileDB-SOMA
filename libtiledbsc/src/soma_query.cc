#include "tiledbsc/soma_query.h"
#include "tiledbsc/common.h"
#include "tiledbsc/logger_public.h"

#include "thread_pool/thread_pool.h"

namespace tiledbsc {
using namespace tiledb;

SOMAQuery::SOMAQuery(SOMA* soma)
    : ctx_(soma->context()) {
    std::vector<ThreadPool::Task> tasks;
    ThreadPool pool{3};

    tasks.emplace_back(pool.execute([&]() {
        mq_obs_ = std::make_unique<ManagedQuery>(soma->open_array("obs"));
        return Status::Ok();
    }));

    tasks.emplace_back(pool.execute([&]() {
        mq_var_ = std::make_unique<ManagedQuery>(soma->open_array("var"));
        return Status::Ok();
    }));
    tasks.emplace_back(pool.execute([&]() {
        mq_x_ = std::make_unique<ManagedQuery>(soma->open_array("X/data"));
        return Status::Ok();
    }));

    pool.wait_all(tasks).ok();
}

std::optional<ColumnBuffers> SOMAQuery::next_results() {
    // Query is complete, return empty results
    if (empty_ || mq_x_->status() == Query::Status::COMPLETE) {
        return std::nullopt;
    }

    // Abort if an invalid column was selected for the obs or var query.
    if (mq_obs_->is_invalid() || mq_var_->is_invalid()) {
        LOG_DEBUG(
            fmt::format("[SOMAQuery] Abort due to invalid selected column."));
        return std::nullopt;
    }

    if (mq_x_->status() == Query::Status::UNINITIALIZED) {
        std::vector<ThreadPool::Task> tasks;
        ThreadPool pool{2};

        // Submit obs and var query tasks in parallel
        tasks.emplace_back(pool.execute([&]() {
            query_and_select(mq_obs_, "obs_id");
            return Status::Ok();
        }));

        tasks.emplace_back(pool.execute([&]() {
            query_and_select(mq_var_, "var_id");
            return Status::Ok();
        }));

        // Block until obs and var tasks complete
        pool.wait_all(tasks).ok();

        // Return empty results if obs or var query was empty
        if (!mq_obs_->total_num_cells() || !mq_var_->total_num_cells()) {
            LOG_DEBUG(fmt::format(
                "Obs or var query was empty: obs num cells={} var num cells={}",
                mq_obs_->total_num_cells(),
                mq_var_->total_num_cells()));
            empty_ = true;
            return std::nullopt;
        }
    }

    // Submit X query
    auto num_cells = mq_x_->submit();
    LOG_DEBUG(fmt::format("*** X cells read = {}", num_cells));

    return mq_x_->results();
}

void SOMAQuery::query_and_select(
    std::unique_ptr<ManagedQuery>& mq, const std::string& dim_name) {
    // Select the dimension column, if all columns are not selected already.
    mq->select_columns({dim_name}, true);

    while (!mq->is_complete()) {
        // Submit the query.
        auto num_cells = mq->submit();
        LOG_DEBUG(fmt::format("*** {} cells read = {}", dim_name, num_cells));

        // If the mq query was sliced and the read returned results, apply the
        // results to the X query.
        // TODO: [optimize] Skip if sliced query returns all results
        if (num_cells && mq->is_sliced()) {
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
    }

    LOG_DEBUG(fmt::format(
        "*** {} total cells read = {}", dim_name, mq->total_num_cells()));
}

}  // namespace tiledbsc
