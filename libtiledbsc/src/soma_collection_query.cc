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
        soma_queries_[name] = soma->query();
    }
}

std::optional<std::unordered_map<std::string, std::shared_ptr<ColumnBuffer>>>
SOMACollectionQuery::next_results() {
    submitted_ = true;

    LOG_DEBUG(fmt::format("[SOMACollectionQuery] Start queries."));

    std::vector<ThreadPool::Task> tasks;
    ThreadPool pool{threads_};

    for (auto& [name, sq] : soma_queries_) {
        LOG_DEBUG(
            fmt::format("[SOMACollectionQuery] Queue query for {}", name));
        tasks.emplace_back(pool.execute([&]() {
            sq->next_results();
            return Status::Ok();
        }));
    }

    // Block until all tasks complete
    pool.wait_all(tasks).ok();

    LOG_DEBUG(fmt::format("[SOMACollectionQuery] Queries done."));

    return std::nullopt;
}

}  // namespace tiledbsc
