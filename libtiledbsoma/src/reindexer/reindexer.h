/**
 * @file   reindexer.h
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
 * This file defines the a pybind11 api into SOMA C++ library.
 */

#ifndef TILEDBSOMA_REINDEXER_H
#define TILEDBSOMA_REINDEXER_H

#include <assert.h>
#include <memory>
#include <stdexcept>
#include <vector>

struct kh_m64_s;

namespace tiledbsoma {

class SOMAContext;

class IntIndexer {
   public:
    /**
     * Perform intitalization of hash and threadpool
     * @param keys pointer to key array of 64bit integers
     * @param size yhr number of keys in the put
     * @param threads number of threads in the thread pool
     */
    void map_locations(const int64_t* keys, size_t size);
    void map_locations(const std::vector<int64_t>& keys) {
        map_locations(keys.data(), keys.size());
    }
    /**
     * Used for parallel lookup using khash
     * @param keys array of keys to lookup
     * @param result array for lookup results
     * @param size // Number of key array
     * @return and array of looked up value (same size as keys)
     */
    void lookup(const int64_t* keys, int64_t* results, size_t size);
    void lookup(
        const std::vector<int64_t>& keys, std::vector<int64_t>& results) {
        if (keys.size() != results.size())
            throw std::runtime_error(
                "The size of input and results arrays must be the same.");

        lookup(keys.data(), results.data(), keys.size());
    }
    IntIndexer(){};
    IntIndexer(std::shared_ptr<tiledbsoma::SOMAContext> context)
        : context_(context) {
    }
    virtual ~IntIndexer();

   private:
    /*
     * The created 64bit hash table
     */
    kh_m64_s* hash_;

    std::shared_ptr<SOMAContext> context_ = nullptr;
    /*
     * Number of elements in the map set by map_locations
     */
    size_t map_size_ = 0;
};

}  // namespace tiledbsoma

#endif  // TILEDBSOMA_REINDEXER_H
