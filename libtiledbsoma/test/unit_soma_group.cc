/**
 * @file   unit_soma_group.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
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
 * This file manages unit tests for the SOMAGroup class
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

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>
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

std::tuple<std::string, uint64_t> create_array(
    const std::string& uri_in,
    Context& ctx,
    int num_cells_per_fragment = 10,
    int num_fragments = 1,
    bool overlap = false,
    bool allow_duplicates = false,
    uint64_t timestamp = 1,
    bool reuse_existing = false) {
    std::string uri = fmt::format(
        "{}-{}-{}-{}-{}",
        uri_in,
        num_cells_per_fragment,
        num_fragments,
        overlap,
        allow_duplicates);
    // Create array, if not reusing the existing array
    if (!reuse_existing) {
        auto vfs = VFS(ctx);
        if (vfs.is_dir(uri)) {
            vfs.remove_dir(uri);
        }

        // Create schema
        ArraySchema schema(ctx, TILEDB_SPARSE);

        auto dim = Dimension::create<int64_t>(
            ctx, "d0", {0, std::numeric_limits<int64_t>::max() - 1});

        Domain domain(ctx);
        domain.add_dimension(dim);
        schema.set_domain(domain);

        auto attr = Attribute::create<int>(ctx, "a0");
        schema.add_attribute(attr);
        schema.set_allows_dups(allow_duplicates);
        schema.check();

        // Create array
        Array::create(uri, schema);
    }

    // Open array for writing
    Array array(ctx, uri, TILEDB_WRITE, timestamp);
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

    if (allow_duplicates) {
        return {uri, nnz};
    }

    // Adjust nnz when overlap is enabled
    if (overlap) {
        nnz = (num_fragments + 1) / 2 * num_cells_per_fragment;
    }

    return {uri, nnz};
}

};  // namespace

TEST_CASE("SOMAGroup: basic") {
    auto ctx = std::make_shared<Context>();

    std::string uri_main_group = "mem://main-group";
    Group::create(*ctx, uri_main_group);

    std::string uri_sub_group = "mem://sub-group";
    Group::create(*ctx, uri_sub_group);

    auto [uri_sub_array, expected_nnz] = create_array("mem://sub-array", *ctx);

    auto soma_group = SOMAGroup::open(
        TILEDB_WRITE, ctx, uri_main_group, "metadata", 1);
    soma_group->add_member(uri_sub_group, false, "subgroup");
    soma_group->add_member(uri_sub_array, false, "subarray");
    soma_group->close();

    soma_group->open(TILEDB_READ, 1);
    REQUIRE(soma_group->ctx() == ctx);
    REQUIRE(soma_group->uri() == uri_main_group);
    REQUIRE(soma_group->get_length() == 2);
    std::map<std::string, std::string> expected_map{
        {"subgroup", uri_sub_group}, {"subarray", uri_sub_array}};
    REQUIRE(expected_map == soma_group->member_to_uri_mapping());
    REQUIRE(soma_group->get_member("subgroup").type() == Object::Type::Group);
    REQUIRE(soma_group->get_member("subarray").type() == Object::Type::Array);
    soma_group->close();

    soma_group->open(TILEDB_WRITE, 3);
    soma_group->remove_member("subgroup");
    soma_group->close();

    soma_group->open(TILEDB_READ, 4);
    REQUIRE(soma_group->get_length() == 1);
    REQUIRE(soma_group->has_member("subgroup") == false);
    REQUIRE(soma_group->has_member("subarray") == true);
    soma_group->close();
}

TEST_CASE("SOMAGroup: metadata") {
    auto ctx = std::make_shared<Context>();

    std::string uri = "mem://unit-test-group";
    Group::create(*ctx, uri);
    auto soma_group = SOMAGroup::open(TILEDB_WRITE, ctx, uri, "metadata", 1);
    int32_t val = 100;
    soma_group->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_group->close();

    soma_group->open(TILEDB_READ, 1);
    REQUIRE(soma_group->has_metadata("md") == true);
    REQUIRE(soma_group->metadata_num() == 1);

    auto mdval = soma_group->get_metadata(0);
    REQUIRE(std::get<MetadataInfo::key>(mdval) == "md");
    REQUIRE(std::get<MetadataInfo::dtype>(mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(mdval)) == 100);

    mdval = soma_group->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::key>(mdval) == "md");
    REQUIRE(std::get<MetadataInfo::dtype>(mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(mdval)) == 100);
    soma_group->close();

    soma_group->open(TILEDB_WRITE, 2);
    soma_group->delete_metadata("md");
    soma_group->close();

    soma_group->open(TILEDB_READ, 3);
    REQUIRE(soma_group->has_metadata("md") == false);
    REQUIRE(soma_group->metadata_num() == 0);
    soma_group->close();
}