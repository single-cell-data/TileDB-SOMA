/**
 * @file   reindexer.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
        std::format("[Re-indexer] Start of Map locations with {} keys", size));
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
    LOG_DEBUG(std::format("[Re-indexer] khash size = {}", hsize));

    LOG_DEBUG(
        std::format("[Re-indexer] Thread pool started and hash table created"));
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
    LOG_DEBUG(std::format(
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
        LOG_DEBUG(std::format(
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
        LOG_DEBUG(std::format(
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
