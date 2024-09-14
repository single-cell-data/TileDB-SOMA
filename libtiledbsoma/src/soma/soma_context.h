/**
 * @file   soma_context.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
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
