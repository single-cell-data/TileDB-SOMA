/**
 * @file   soma_collection_query.cc
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
 * This file defines the query class for a soma collection object.
 */

#include "tiledbsc/soma_collection_query.h"
#include "tiledbsc/logger_public.h"
#include "tiledbsc/soma_collection.h"

#include "thread_pool/thread_pool.h"

namespace tiledbsc {
using namespace tiledb;

SOMACollectionQuery::SOMACollectionQuery(SOMACollection* soco) {
    std::vector<ThreadPool::Task> tasks;
    ThreadPool pool{threads_};

    for (auto& [name, soma] : soco->get_somas()) {
        LOG_DEBUG(fmt::format("Get SOMA query: {}", name));
        soma_queries_.push_back(soma->query(name));
    }
}

std::optional<MultiArrayBuffers> SOMACollectionQuery::next_results() {
    submitted_ = true;

    LOG_DEBUG(fmt::format("[SOMACollectionQuery] Start queries."));

    std::vector<ThreadPool::Task> tasks;
    ThreadPool pool{threads_};

    for (auto& sq : soma_queries_) {
        if (!sq->is_complete()) {
            LOG_DEBUG(fmt::format(
                "[SOMACollectionQuery] Queue query for {}", sq->name()));
            tasks.emplace_back(pool.execute([&]() {
                sq->next_results();
                return Status::Ok();
            }));
        }
    }

    // Block until all tasks complete
    pool.wait_all(tasks).ok();

    LOG_DEBUG(fmt::format("[SOMACollectionQuery] Queries done."));

    MultiArrayBuffers results;

    for (auto& sq : soma_queries_) {
        if (sq->results().has_value()) {
            LOG_DEBUG(fmt::format(
                "[SOMACollectionQuery] SOMA {} has {} results.",
                sq->name(),
                sq->results()->begin()->second.begin()->second->size()));

            // Merge soma query results into results to be returned
            results.merge(*sq->results());
        }
    }

    if (results.empty()) {
        return std::nullopt;
    }
    return results;
}

}  // namespace tiledbsc
