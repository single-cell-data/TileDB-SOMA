/**
 * @file   soma_context.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAContext class.
 */
#include <thread>

#include <thread_pool/thread_pool.h>
#include "../utils/common.h"
#include "../utils/logger.h"
#include "soma_context.h"

namespace tiledbsoma {

std::shared_ptr<ThreadPool>& SOMAContext::thread_pool() {
    const std::lock_guard<std::mutex> lock(thread_pool_mutex_);
    // The first thread that gets here will create the context thread pool
    if (thread_pool_ == nullptr) {
        auto cfg = tiledb_config();
        auto concurrency = std::thread::hardware_concurrency();
        if (cfg.find(CONFIG_KEY_COMPUTE_CONCURRENCY_LEVEL) != cfg.end()) {
            auto value_str = cfg[CONFIG_KEY_COMPUTE_CONCURRENCY_LEVEL];
            try {
                concurrency = std::stoull(value_str);
            } catch (const std::exception& e) {
                throw TileDBSOMAError(std::format(
                    "[SOMAContext] Error parsing {}: '{}' ({}) - must be a "
                    "postive integer.",
                    CONFIG_KEY_COMPUTE_CONCURRENCY_LEVEL,
                    value_str,
                    e.what()));
            }
        }

        int thread_count = std::min(std::max(1u, concurrency), 1024u);
        thread_pool_ = std::make_shared<ThreadPool>(thread_count);
    }
    return thread_pool_;
}
}  // namespace tiledbsoma
