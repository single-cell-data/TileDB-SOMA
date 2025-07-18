/**
 * @file   test_indexer.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for re-indexer
 */

#include <reindexer/reindexer.h>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <string>
#include <tiledb/tiledb>
#include <unordered_map>
#include <vector>
#include "external/khash/khash.h"

namespace {

/**
 * @brief Create multiple 64 bit vector and try the reindexer compared to khash
 * as baseline.
 *
 */

KHASH_MAP_INIT_INT64(m64, int64_t)

// Identify the test in which keys are not unique and we should expect
// exception.
std::vector<bool> uniqueness = {false, false, true, true, true};

bool run_test(size_t id, std::vector<int64_t> keys, std::vector<int64_t> lookups) {
    try {
        std::vector<int64_t> indexer_results;
        indexer_results.resize(lookups.size());

        tiledbsoma::IntIndexer indexer;
        indexer.map_locations(keys);
        auto* hash = kh_init(m64);
        int ret;
        khint64_t k;

        for (size_t i = 0; i < keys.size(); i++) {
            k = kh_put(m64, hash, keys[i], &ret);
            assert(k != kh_end(hash));
            kh_val(hash, k) = i;
        }

        indexer.lookup(lookups, indexer_results);
        std::vector<int64_t> kh_results;
        kh_results.resize(lookups.size());
        for (size_t i = 0; i < lookups.size(); i++) {
            auto k = kh_get(m64, hash, lookups[i]);
            if (k == kh_end(hash)) {
                // According to pandas behavior
                kh_results[i] = -1;
            } else {
                kh_results[i] = kh_val(hash, k);
            }
        }
        kh_destroy(m64, hash);
        return indexer_results == kh_results;
    } catch (const std::runtime_error& e) {
        // Path with non unique keys
        if (!uniqueness[id]) {
            return true;
        }
        return false;
    }
}

// Test data
std::vector<std::unordered_map<std::string, std::vector<int64_t>>> test_data = {
    {
        {"keys", {-1, -1, -1, 0, 0, 0}},
        {"lookups", {1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5}},
    },
    {{"keys",
      {
          -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5,
      }},
     {"lookups", {-10000, 1, 2, 3, 5, 6}}},
    {
        {"keys", {-1, 1, 2, 3, 4, 5}},
        {"lookups",
         {
             -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5,
         }},
    },
    {
        {"keys", {-10000, -100000, 200000, 5, 1, 7}},
        {"lookups",
         {
             -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5,
         }},
    },
    {
        {"keys", {-10000, -200000, 1000, 3000, 1, 2}},
        {"lookups",
         {
             -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5, -1, 1, 2, 3, 4, 5,
         }},
    }};

TEST_CASE("C++ re-indexer") {
    for (size_t test = 0; test < test_data.size(); test++) {
        bool result = run_test(test, test_data[test]["keys"], test_data[test]["lookups"]);
        if (!result) {
            throw std::runtime_error("Test " + std::to_string(test) + " failed");
        }
    }
}
}  // namespace
