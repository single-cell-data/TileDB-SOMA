/**
 * @file   unit_managed_query.cc
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
 * This file manages unit tests for managed queries
 */

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_predicate.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
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

bool VERBOSE = false;

auto create_array(const std::string& uri, Context& ctx) {
    // Delete array if it exists
    auto vfs = VFS(ctx);
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    // Create schema
    ArraySchema schema(ctx, TILEDB_SPARSE);
    auto dim = Dimension::create(
        ctx, "d0", TILEDB_STRING_ASCII, nullptr, nullptr);
    dim.set_cell_val_num(TILEDB_VAR_NUM);

    Domain domain(ctx);
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<std::string>(ctx, "a0");
    attr.set_nullable(true);
    schema.add_attribute(attr);
    schema.check();

    // Create array and open for writing
    Array::create(uri, schema);
    Array array(ctx, uri, TILEDB_WRITE);

    std::vector<std::string> d0 = {
        "a", "bb", "ccc", "dddd", "eeeee", "fffffff"};
    auto [d0_data, d0_offsets] = util::to_varlen_buffers(d0, false);

    std::vector<std::string> a0 = {
        "A", "BB", "CCC", "DDDD", "EEEEE", "FFFFFFF"};
    auto [a0_data, a0_offsets] = util::to_varlen_buffers(a0, false);

    std::vector<uint8_t> a0_valids = {1, 1, 1, 1, 0, 0};

    // Write data to array and close the array
    Query query(ctx, array);
    query.set_layout(TILEDB_UNORDERED)
        .set_data_buffer("d0", d0_data)
        .set_offsets_buffer("d0", d0_offsets)
        .set_data_buffer("a0", a0_data)
        .set_offsets_buffer("a0", a0_offsets)
        .set_validity_buffer("a0", a0_valids);
    query.submit();
    array.close();

    // Open the array for reading and return a shared pointer
    return std::tuple(
        std::make_shared<Array>(ctx, uri, TILEDB_READ), d0, a0, a0_valids);
}

};  // namespace

TEST_CASE("ManagedQuery: Basic execution test") {
    if (VERBOSE) {
        LOG_CONFIG("debug");
    }

    std::string uri = "mem://unit-test-array";
    auto ctx = Context();
    auto [array, d0, a0, _] = create_array(uri, ctx);

    auto mq = ManagedQuery(array);
    mq.submit();

    auto results = mq.results();
    REQUIRE(mq.results_complete());

    auto num_cells = mq.total_num_cells();
    REQUIRE(num_cells == d0.size());

    REQUIRE_THAT(d0, Equals(mq.strings("d0")));
    REQUIRE_THAT(a0, Equals(mq.strings("a0")));
}

TEST_CASE("ManagedQuery: Select test") {
    std::string uri = "mem://unit-test-array";
    auto ctx = Context();
    auto [array, d0, a0, _] = create_array(uri, ctx);

    auto mq = ManagedQuery(array);
    mq.select_columns({"a0"});
    mq.select_points<std::string>("d0", {"a"});
    mq.submit();

    auto results = mq.results();
    REQUIRE(mq.results_complete());

    auto num_cells = mq.total_num_cells();
    REQUIRE(num_cells == 1);

    REQUIRE_THROWS(mq.data<int>("a1"));
    REQUIRE_THROWS(mq.strings("d0"));
    REQUIRE_THROWS(mq.string_view("d1", 0));

    REQUIRE_THAT(
        std::string(a0[0]), Equals(std::string(mq.string_view("a0", 0))));
}

TEST_CASE("ManagedQuery: Validity test") {
    std::string uri = "mem://unit-test-array";
    auto ctx = Context();
    auto [array, d0, a0, a0_valids] = create_array(uri, ctx);

    auto mq = ManagedQuery(array);
    mq.submit();

    auto results = mq.results();
    REQUIRE(mq.results_complete());

    auto num_cells = mq.total_num_cells();
    REQUIRE(num_cells == d0.size());

    // Convert span to vector
    auto valids = mq.validity("a0");
    std::vector<uint8_t> a0_valids_actual;
    a0_valids_actual.assign(valids.begin(), valids.end());

    REQUIRE_THAT(d0, Equals(mq.strings("d0")));
    REQUIRE_THAT(a0, Equals(mq.strings("a0")));
    REQUIRE_THAT(a0_valids, Equals(a0_valids_actual));
}
