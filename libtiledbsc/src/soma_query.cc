/**
 * @file   soma_query.cc
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
 * This file defines the query class for a soma object.
 */

#include "tiledbsc/soma_query.h"
#include "tiledbsc/common.h"
#include "tiledbsc/logger_public.h"

#include "thread_pool/thread_pool.h"

namespace tiledbsc {
using namespace tiledb;

SOMAQuery::SOMAQuery(SOMA* soma, std::string name)
    : name_(name)
    , ctx_(soma->context()) {
    std::vector<ThreadPool::Task> tasks;
    ThreadPool pool{3};

    tasks.emplace_back(pool.execute([&]() {
        std::string name{"obs"};
        std::string unique_name{name_};
        if (!name_.empty()) {
            unique_name += "/" + name;
        }
        mq_obs_ = std::make_unique<ManagedQuery>(
            soma->open_array(name), unique_name);
        return Status::Ok();
    }));

    tasks.emplace_back(pool.execute([&]() {
        std::string name{"var"};
        std::string unique_name{name_};
        if (!name_.empty()) {
            unique_name += "/" + name;
        }
        mq_var_ = std::make_unique<ManagedQuery>(
            soma->open_array(name), unique_name);
        return Status::Ok();
    }));

    tasks.emplace_back(pool.execute([&]() {
        std::string name{"X/data"};
        std::string unique_name{name_};
        if (!name_.empty()) {
            unique_name += "/" + name;
        }
        mq_x_ = std::make_unique<ManagedQuery>(
            soma->open_array(name), unique_name);
        return Status::Ok();
    }));

    pool.wait_all(tasks).ok();
}

std::optional<MultiArrayBuffers> SOMAQuery::next_results() {
    // Query is complete, return empty results
    if (empty_ || mq_x_->status() == Query::Status::COMPLETE) {
        results_.clear();
        complete_ = true;
        return std::nullopt;
    }

    // Abort if an invalid column was selected for the obs or var query.
    if (mq_obs_->is_invalid() || mq_var_->is_invalid()) {
        LOG_DEBUG(fmt::format(
            "[SOMAQuery] [{}] Abort due to invalid selected column.", name_));
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
                "[SOMAQuery] [{}] Obs or var query was empty: obs num cells={} "
                "var num cells={}",
                name_,
                mq_obs_->total_num_cells(),
                mq_var_->total_num_cells()));
            empty_ = true;
            return std::nullopt;
        }
    }

    // Submit X query
    auto num_cells = mq_x_->submit();
    LOG_DEBUG(fmt::format(
        "[SOMAQuery] [{}/X/data] cells read = {}", name_, num_cells));

    // Save results in MultiArrayBuffers
    results_.clear();

    if (!mq_obs_->results().empty()) {
        results_[name_ + "/obs"] = mq_obs_->results();
    }
    if (!mq_var_->results().empty()) {
        results_[name_ + "/var"] = mq_var_->results();
    }
    if (!mq_x_->results().empty()) {
        results_[name_ + "/X/data"] = mq_x_->results();
    }

    return results_;
}

void SOMAQuery::query_and_select(
    std::unique_ptr<ManagedQuery>& mq, const std::string& dim_name) {
    // Select the dimension column, if all columns are not selected already.
    mq->select_columns({dim_name}, true);

    while (!mq->is_complete()) {
        // Submit the query.
        auto num_cells = mq->submit();
        LOG_DEBUG(fmt::format(
            "[SOMAQuery] [{}] {} cells read = {}", name_, dim_name, num_cells));

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
        "[SOMAQuery] [{}] {} total cells read = {}",
        name_,
        dim_name,
        mq->total_num_cells()));
}

}  // namespace tiledbsc
