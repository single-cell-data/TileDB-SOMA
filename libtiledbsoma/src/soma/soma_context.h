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
 * This file defines the SOMAContext class.
 */

#ifndef SOMA_CONTEXT
#define SOMA_CONTEXT

#include <map>
#include <mutex>
#include <string>
#include <tiledb/tiledb>

namespace tiledbsoma {
class ThreadPool;

using namespace tiledb;

class SOMAContext {
    // Controls concurrency level for SOMA compute thread pool. Defaults to host
    // CPU count.
    inline static const std::string
        CONFIG_KEY_COMPUTE_CONCURRENCY_LEVEL = "soma.compute_concurrency_level";

   public:
    //===================================================================
    //= public non-static
    //===================================================================
    SOMAContext()
        : ctx_(std::make_shared<Context>(Config({})))
        , thread_pool_mutex_(){};

    SOMAContext(std::map<std::string, std::string> tiledb_config)
        : ctx_(std::make_shared<Context>(Config(tiledb_config)))
        , thread_pool_mutex_(){};

    bool operator==(const SOMAContext& other) const {
        return ctx_ == other.ctx_;
    }

    std::shared_ptr<Context> tiledb_ctx() const {
        return ctx_;
    }

    std::map<std::string, std::string> tiledb_config() const {
        std::map<std::string, std::string> cfg;
        for (auto& it : ctx_->config())
            cfg[it.first] = it.second;
        return cfg;
    }

    std::shared_ptr<ThreadPool>& thread_pool();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // Threadpool
    std::shared_ptr<ThreadPool> thread_pool_ = nullptr;

    // Semaphore to create and use the thread_pool
    std::mutex thread_pool_mutex_;
};
}  // namespace tiledbsoma

#endif  // SOMA_CONTEXT
