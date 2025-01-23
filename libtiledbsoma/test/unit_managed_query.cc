/**
 * @file   unit_managed_query.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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

#include <tiledb/tiledb>
#include <tiledbsoma/tiledbsoma>
#include "utils/util.h"

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

    std::string dim_name = "d0";
    std::string attr_name = "a0";

    // Create schema
    ArraySchema schema(ctx, TILEDB_SPARSE);
    auto dim = Dimension::create(
        ctx, dim_name, TILEDB_STRING_ASCII, nullptr, nullptr);
    dim.set_cell_val_num(TILEDB_VAR_NUM);

    Domain domain(ctx);
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<std::string>(ctx, attr_name);
    attr.set_nullable(true);
    schema.add_attribute(attr);
    schema.check();

    // Create array and open for writing
    Array::create(uri, std::move(schema));
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
        .set_data_buffer(dim_name, d0_data)
        .set_offsets_buffer(dim_name, d0_offsets)
        .set_data_buffer(attr_name, a0_data)
        .set_offsets_buffer(attr_name, a0_offsets)
        .set_validity_buffer(attr_name, a0_valids);
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
    std::string dim_name = "d0";
    std::string attr_name = "a0";

    auto ctx = std::make_shared<Context>();
    auto [array, d0, a0, _] = create_array(uri, *ctx);

    auto mq = ManagedQuery(array, ctx);
    auto results = mq.read_next();
    REQUIRE(mq.results_complete());

    auto num_cells = mq.total_num_cells();
    REQUIRE(num_cells == d0.size());

    REQUIRE_THAT(d0, Equals(mq.strings(dim_name)));
    REQUIRE_THAT(a0, Equals(mq.strings(attr_name)));
}

TEST_CASE("ManagedQuery: Select test") {
    std::string uri = "mem://unit-test-array";
    std::string dim_name = "d0";
    std::string attr_name = "a0";

    auto ctx = std::make_shared<Context>();
    auto [array, d0, a0, _] = create_array(uri, *ctx);

    auto mq = ManagedQuery(array, ctx);
    mq.select_columns({attr_name});
    mq.select_points<std::string>(dim_name, {"a"});
    auto results = mq.read_next();
    REQUIRE(mq.results_complete());

    auto num_cells = mq.total_num_cells();
    REQUIRE(num_cells == 1);

    REQUIRE_THROWS(mq.data<int32_t>("a1"));
    REQUIRE_THROWS(mq.strings(dim_name));
    REQUIRE_THROWS(mq.string_view("d1", 0));

    REQUIRE_THAT(
        std::string(a0[0]), Equals(std::string(mq.string_view(attr_name, 0))));
}

TEST_CASE("ManagedQuery: Validity test") {
    std::string uri = "mem://unit-test-array";
    std::string dim_name = "d0";
    std::string attr_name = "a0";

    auto ctx = std::make_shared<Context>();
    auto [array, d0, a0, a0_valids] = create_array(uri, *ctx);

    auto mq = ManagedQuery(array, ctx);
    auto results = mq.read_next();
    REQUIRE(mq.results_complete());

    auto num_cells = mq.total_num_cells();
    REQUIRE(num_cells == d0.size());

    // Convert span to vector
    auto valids = mq.validity(attr_name);
    std::vector<uint8_t> a0_valids_actual;
    a0_valids_actual.assign(valids.begin(), valids.end());

    REQUIRE_THAT(d0, Equals(mq.strings(dim_name)));
    REQUIRE_THAT(a0, Equals(mq.strings(attr_name)));
    REQUIRE_THAT(a0_valids, Equals(a0_valids_actual));
}
