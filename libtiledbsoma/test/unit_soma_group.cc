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
    const std::string& uri,
    Context& ctx,
    int num_cells_per_fragment = 10,
    int num_fragments = 1,
    bool overlap = false,
    bool allow_duplicates = false,
    uint64_t timestamp = 1,
    bool reuse_existing = false) {
    std::string dim_name = "d0";
    std::string attr_name = "a0";

    // Create array, if not reusing the existing array
    if (!reuse_existing) {
        auto vfs = VFS(ctx);
        if (vfs.is_dir(uri)) {
            vfs.remove_dir(uri);
        }

        // Create schema
        ArraySchema schema(ctx, TILEDB_SPARSE);

        auto dim = Dimension::create<int64_t>(
            ctx, dim_name, {0, std::numeric_limits<int64_t>::max() - 1});

        Domain domain(ctx);
        domain.add_dimension(dim);
        schema.set_domain(domain);

        auto attr = Attribute::create<int32_t>(ctx, attr_name);
        schema.add_attribute(attr);
        schema.set_allows_dups(allow_duplicates);
        schema.check();

        // Create array
        Array::create(uri, std::move(schema));
    }

    // Open array for writing
    Array array(ctx, uri, TILEDB_WRITE, TemporalPolicy(TimeTravel, timestamp));
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
            .set_data_buffer(dim_name, d0)
            .set_data_buffer(attr_name, a0);
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
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri_main_group = "mem://main-group";
    SOMAGroup::create(ctx, uri_main_group, "NONE");

    std::string uri_sub_group = "mem://sub-group";
    SOMAGroup::create(ctx, uri_sub_group, "NONE");

    auto [uri_sub_array, expected_nnz] = create_array(
        "mem://sub-array", *ctx->tiledb_ctx());

    auto soma_group = SOMAGroup::open(
        OpenMode::write, uri_main_group, ctx, "metadata", TimestampRange(0, 1));
    soma_group->set(uri_sub_group, URIType::absolute, "subgroup", "SOMAGroup");
    soma_group->set(uri_sub_array, URIType::absolute, "subarray", "SOMAArray");
    soma_group->close();

    std::map<std::string, SOMAGroupEntry> expected_map{
        {"subgroup", SOMAGroupEntry(uri_sub_group, "SOMAGroup")},
        {"subarray", SOMAGroupEntry(uri_sub_array, "SOMAArray")}};

    soma_group->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(soma_group->ctx() == ctx);
    REQUIRE(soma_group->uri() == uri_main_group);
    REQUIRE(soma_group->count() == 2);
    REQUIRE(expected_map == soma_group->members_map());
    REQUIRE(soma_group->get("subgroup").type() == Object::Type::Group);
    REQUIRE(soma_group->get("subarray").type() == Object::Type::Array);
    soma_group->close();

    soma_group->open(OpenMode::write, TimestampRange(0, 3));
    REQUIRE(expected_map == soma_group->members_map());
    soma_group->del("subgroup");
    soma_group->close();

    soma_group->open(OpenMode::read, TimestampRange(0, 4));
    REQUIRE(soma_group->count() == 1);
    REQUIRE(soma_group->has("subgroup") == false);
    REQUIRE(soma_group->has("subarray") == true);
    soma_group->close();
}

TEST_CASE("SOMAGroup: metadata") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-group";
    SOMAGroup::create(ctx, uri, "NONE", TimestampRange(0, 2));
    auto soma_group = SOMAGroup::open(
        OpenMode::write, uri, ctx, "metadata", TimestampRange(1, 1));
    int32_t val = 100;
    soma_group->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_group->close();

    // Read metadata
    soma_group->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(soma_group->metadata_num() == 3);
    REQUIRE(soma_group->has_metadata("soma_object_type"));
    REQUIRE(soma_group->has_metadata("soma_encoding_version"));
    REQUIRE(soma_group->has_metadata("md"));
    auto mdval = soma_group->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_group->close();

    // md should not be available at (2, 2)
    soma_group->open(OpenMode::read, TimestampRange(2, 2));
    REQUIRE(soma_group->metadata_num() == 2);
    REQUIRE(soma_group->has_metadata("soma_object_type"));
    REQUIRE(soma_group->has_metadata("soma_encoding_version"));
    REQUIRE(!soma_group->has_metadata("md"));
    soma_group->close();

    // Metadata should also be retrievable in write mode
    soma_group->open(OpenMode::write, TimestampRange(0, 2));
    REQUIRE(soma_group->metadata_num() == 3);
    REQUIRE(soma_group->has_metadata("soma_object_type"));
    REQUIRE(soma_group->has_metadata("soma_encoding_version"));
    REQUIRE(soma_group->has_metadata("md"));
    mdval = soma_group->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    soma_group->delete_metadata("md");
    mdval = soma_group->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_group->close();

    // Confirm delete in read mode
    soma_group->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(!soma_group->has_metadata("md"));
    REQUIRE(soma_group->metadata_num() == 2);
}

TEST_CASE("SOMAGroup: dataset_type") {
    auto ctx = std::make_shared<SOMAContext>();
    SOMAGroup::create(ctx, "mem://experiment", "SOMAExperiment");
    SOMAGroup::create(ctx, "mem://collection", "SOMACollection");
    SOMAGroup::create(ctx, "mem://measurement", "SOMAMeasurement");

    auto experiment = SOMAGroup::open(OpenMode::read, "mem://experiment", ctx);
    auto collection = SOMAGroup::open(OpenMode::read, "mem://collection", ctx);
    auto measurement = SOMAGroup::open(
        OpenMode::read, "mem://measurement", ctx);

    REQUIRE(!collection->has_metadata("dataset_type"));
    REQUIRE(!measurement->has_metadata("dataset_type"));

    REQUIRE(experiment->has_metadata("dataset_type"));

    std::string expect = "soma";

    // tuple of dtype, count, void*:
    auto dataset_type = experiment->get_metadata("dataset_type");
    auto bytes = (const char*)std::get<MetadataInfo::value>(*dataset_type);
    auto count = std::get<MetadataInfo::num>(*dataset_type);
    auto actual = std::string(bytes, count);

    REQUIRE(actual == expect);
}
