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
