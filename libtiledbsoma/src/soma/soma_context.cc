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
 *   This file defines the SOMAContext class.
 */
#include "soma_context.h"
#include <thread_pool/thread_pool.h>

namespace tiledbsoma {

std::shared_ptr<ThreadPool>& SOMAContext::thread_pool() {
    const std::lock_guard<std::mutex> lock(thread_pool_mutex_);
    // The first thread that gets here will create the context thread pool
    if (thread_pool_ == nullptr) {
        auto cfg = tiledb_config();
        int concurrency = 10;
        if (cfg.find("sm.compute_concurrency_level") != cfg.end()) {
            concurrency = std::stoi(cfg["sm.compute_concurrency_level"]);
        }
        int thread_count = std::max(1, concurrency / 2);
        if (thread_count > 1) {
            thread_pool_ = std::make_shared<ThreadPool>(thread_count);
        }
    }
    return thread_pool_;
}
}  // namespace tiledbsoma
