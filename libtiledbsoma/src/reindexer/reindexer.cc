/**
 * @file   reindexer.cc
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
 * This file defines the a pybind11 api in the SOMA C++ library.
 */

#include "reindexer.h"
#include <thread_pool/thread_pool.h>
#include <thread>
#include "khash.h"
#include "soma/enums.h"
#include "soma/soma_context.h"
#include "utils/arrow_adapter.h"
#include "utils/common.h"
#include "utils/logger.h"

// Typedef for a 64-bit khash table
KHASH_MAP_INIT_INT64(m64, int64_t)

namespace tiledbsoma {

void IntIndexer::map_locations(const int64_t* keys, size_t size) {
    map_size_ = size;

    // Handling edge cases
    if (size == 0) {
        return;
    }

    hash_ = kh_init(m64);
    kh_resize(m64, hash_, size * 1.25);
    int ret;
    khint64_t k;
    int64_t counter = 0;
    // Hash map construction
    LOG_DEBUG(
        fmt::format("[Re-indexer] Start of Map locations with {} keys", size));
    for (size_t i = 0; i < size; i++) {
        k = kh_put(m64, hash_, keys[i], &ret);
        assert(k != kh_end(hash_));
        kh_val(hash_, k) = counter;
        counter++;
    }
    if (kh_size(hash_) != size) {
        throw std::runtime_error("There are duplicate keys.");
    }
    auto hsize = kh_size(hash_);
    LOG_DEBUG(fmt::format("[Re-indexer] khash size = {}", hsize));

    LOG_DEBUG(
        fmt::format("[Re-indexer] Thread pool started and hash table created"));
}

void IntIndexer::lookup(const int64_t* keys, int64_t* results, size_t size) {
    if (size == 0) {
        return;
    }
    // Single thread checks
    if (context_ == nullptr || context_->thread_pool() == nullptr ||
        context_->thread_pool()->concurrency_level() == 1) {
        for (size_t i = 0; i < size; i++) {
            auto k = kh_get(m64, hash_, keys[i]);
            if (k == kh_end(hash_)) {
                // According to pandas behavior
                results[i] = -1;
            } else {
                results[i] = kh_val(hash_, k);
            }
        }
        return;
    }
    LOG_DEBUG(fmt::format(
        "Lookup with thread concurrency {} on data size {}",
        context_->thread_pool()->concurrency_level(),
        size));

    std::vector<tiledbsoma::ThreadPool::Task> tasks;

    size_t thread_chunk_size = size /
                               context_->thread_pool()->concurrency_level();
    if (thread_chunk_size == 0) {
        thread_chunk_size = 1;
    }

    for (size_t i = 0; i < size; i += thread_chunk_size) {
        size_t start = i;
        size_t end = i + thread_chunk_size;
        if (end > size) {
            end = size;
        }
        LOG_DEBUG(fmt::format(
            "Creating tileDB task for the range from {} to {} ", start, end));
        tiledbsoma::ThreadPool::Task task = context_->thread_pool()->execute(
            [this, start, end, &results, &keys]() {
                for (size_t i = start; i < end; i++) {
                    auto k = kh_get(m64, hash_, keys[i]);
                    if (k == kh_end(hash_)) {
                        // According to pandas behavior
                        results[i] = -1;
                    } else {
                        results[i] = kh_val(hash_, k);
                    }
                }
                return tiledbsoma::Status::Ok();
            });
        assert(task.valid());
        tasks.emplace_back(std::move(task));
        LOG_DEBUG(fmt::format(
            "Task for the range from {} to {} inserted in the queue",
            start,
            end));
    }
    context_->thread_pool()->wait_all(tasks);
}

IntIndexer::~IntIndexer() {
    if (map_size_ > 0) {
        kh_destroy(m64, this->hash_);
    }
}

}  // namespace tiledbsoma
