/**
 * @file   indexer.cc
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
 * This file defines the Reindexer bindings.
 */

#include <tiledbsoma/reindexer/reindexer.h>

#include "common.h"

#define DENUM(x) .value(#x, TILEDB_##x)
namespace libtiledbsomacpp {

namespace py = pybind11;
using namespace py::literals;
using namespace tiledbsoma;

void load_reindexer(py::module& m) {
    // Efficient C++ re-indexing (aka hashing unique key values to an index
    // between 0 and number of keys - 1) based on khash
    py::class_<IntIndexer>(m, "IntIndexer")
        .def(py::init<>())
        .def(py::init<std::vector<int64_t>&, int>())
        .def(
            "map_locations",
            [](IntIndexer& indexer,
               py::array_t<int64_t> keys,
               int num_threads) {
                auto buffer = keys.request();
                int64_t* data = static_cast<int64_t*>(buffer.ptr);
                size_t length = buffer.shape[0];
                indexer.map_locations(keys.data(), keys.size(), num_threads);
            })
        .def(
            "map_locations",
            [](IntIndexer& indexer,
               std::vector<int64_t> keys,
               int num_threads) {
                indexer.map_locations(keys.data(), keys.size(), num_threads);
            })
        // Perform lookup for a large input array of keys and return the looked
        // up value array (passing ownership from C++ to python)
        .def(
            "get_indexer",
            [](IntIndexer& indexer, py::array_t<int64_t> lookups) {
                auto input_buffer = lookups.request();
                int64_t* input_ptr = static_cast<int64_t*>(input_buffer.ptr);
                size_t size = input_buffer.shape[0];
                auto results = py::array_t<int64_t>(size);
                auto results_buffer = results.request();
                size_t results_size = results_buffer.shape[0];

                int64_t* results_ptr = static_cast<int64_t*>(
                    results_buffer.ptr);

                indexer.lookup(input_ptr, results_ptr, size);
                return results;
            })
        // Perform lookup for a large input array of keys and writes the looked
        // up values into previously allocated array (works for the cases in
        // which python and R pre-allocate the array)
        .def(
            "get_indexer",
            [](IntIndexer& indexer,
               py::array_t<int64_t> lookups,
               py::array_t<int64_t>& results) {
                auto input_buffer = lookups.request();
                int64_t* input_ptr = static_cast<int64_t*>(input_buffer.ptr);
                size_t size = input_buffer.shape[0];

                auto results_buffer = results.request();
                int64_t* results_ptr = static_cast<int64_t*>(
                    results_buffer.ptr);
                size_t results_size = input_buffer.shape[0];
                indexer.lookup(input_ptr, input_ptr, size);
            });
}
}  // namespace libtiledbsomacpp