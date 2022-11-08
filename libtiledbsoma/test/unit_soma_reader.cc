/**
 * @file   unit_soma_reader.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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
 * This file manages unit tests for the SOMAReader class
 */

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_predicate.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <numeric>
#include <random>

#include <tiledbsoma/util.h>
#include <tiledb/tiledb>
#include <tiledbsoma/tiledbsoma>

using namespace tiledb;
using namespace tiledbsoma;
using namespace Catch::Matchers;

#ifndef TILEDBSOMA_SOURCE_ROOT
#define TILEDBSOMA_SOURCE_ROOT "not_defined"
#endif

const std::string src_path = TILEDBSOMA_SOURCE_ROOT;

namespace {

uint64_t create_array(
    const std::string& uri,
    Context& ctx,
    int num_cells_per_fragment = 10,
    int num_fragments = 1,
    bool overlap = false,
    bool allow_duplicates = false) {
    // Delete array if it exists
    auto vfs = VFS(ctx);
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    // Create schema
    ArraySchema schema(ctx, TILEDB_SPARSE);

    auto dim = Dimension::create<int64_t>(
        ctx, "d0", {0, std::numeric_limits<int64_t>::max()});

    Domain domain(ctx);
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int>(ctx, "a0");
    schema.add_attribute(attr);
    schema.set_allows_dups(allow_duplicates);
    schema.check();

    // Create array and open for writing
    Array::create(uri, schema);
    Array array(ctx, uri, TILEDB_WRITE);
    if (LOG_DEBUG_ENABLED()) {
        array.schema().dump();
    }

    // Generate fragments in random order
    std::vector<int> frags(num_fragments);
    std::iota(frags.begin(), frags.end(), 0);
    std::shuffle(frags.begin(), frags.end(), std::random_device{});

    for (auto i : frags) {
        std::vector<int64_t> d0(num_cells_per_fragment);
        for (int j = 0; j < num_cells_per_fragment; j++) {
            // Overlap odd fragments when generating overlaps
            if (overlap && i & 1) {
                d0[j] = j + num_cells_per_fragment * (i - 1);
            } else {
                d0[j] = j + num_cells_per_fragment * i;
            }
        }
        std::vector<int> a0(num_cells_per_fragment, i);

        // Write data to array
        Query query(ctx, array);
        query.set_layout(TILEDB_UNORDERED)
            .set_data_buffer("d0", d0)
            .set_data_buffer("a0", a0);
        query.submit();
    }

    array.close();

    uint64_t nnz = num_fragments * num_cells_per_fragment;

    // Adjust nnz when overlap is enabled
    if (overlap) {
        nnz = (num_fragments + 1) / 2 * num_cells_per_fragment;
    }

    return nnz;
}

};  // namespace

TEST_CASE("SOMAReader: nnz") {
    auto num_fragments = GENERATE(1, 10);
    auto overlap = GENERATE(false, true);
    auto allow_duplicates = GENERATE(false, true);
    int num_cells_per_fragment = 128;

    SECTION(fmt::format(
        " - fragments={}, overlap={}, allow_duplicates={}",
        num_fragments,
        overlap,
        allow_duplicates)) {
        // Create array
        std::string uri = "mem://unit-test-array";
        auto ctx = std::make_shared<Context>();
        auto expected_nnz = create_array(
            uri,
            *ctx,
            num_cells_per_fragment,
            num_fragments,
            overlap,
            allow_duplicates);

        // Get total cell num
        auto sr = SOMAReader::open(ctx, uri);

        uint64_t nnz;
        if (num_fragments > 1 && overlap && allow_duplicates) {
            REQUIRE_THROWS(nnz = sr->nnz());
        } else {
            nnz = sr->nnz();
            REQUIRE(nnz == expected_nnz);
        }
    }
}
