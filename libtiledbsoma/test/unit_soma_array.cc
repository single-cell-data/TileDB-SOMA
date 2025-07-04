/**
 * @file   unit_soma_array.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMAArray class
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
#include <tiledbsoma/tiledbsoma>
#include "common.h"
#include "utils/util.h"

using namespace tiledb;
using namespace tiledbsoma;
using namespace Catch::Matchers;

#ifndef TILEDBSOMA_SOURCE_ROOT
#define TILEDBSOMA_SOURCE_ROOT "not_defined"
#endif

namespace {

std::tuple<std::string, uint64_t> create_array(
    const std::string& uri,
    std::shared_ptr<SOMAContext> ctx,
    int num_cells_per_fragment = 10,
    int num_fragments = 1,
    bool overlap = false,
    bool allow_duplicates = false) {
    auto vfs = VFS(*ctx->tiledb_ctx());
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    const char* dim_name = "d0";
    const char* attr_name = "a0";

    // Create schema
    ArraySchema schema(*ctx->tiledb_ctx(), TILEDB_SPARSE);

    auto dim = Dimension::create<int64_t>(*ctx->tiledb_ctx(), dim_name, {0, std::numeric_limits<int64_t>::max() - 1});

    Domain domain(*ctx->tiledb_ctx());
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int32_t>(*ctx->tiledb_ctx(), attr_name);
    schema.add_attribute(attr);
    schema.set_allows_dups(allow_duplicates);
    schema.check();

    // Create array
    SOMAArray::create(ctx, uri, std::move(schema), "NONE", std::nullopt, TimestampRange(0, 2));

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

std::tuple<std::vector<int64_t>, std::vector<int32_t>> write_array(
    const std::string& uri,
    std::shared_ptr<SOMAContext> ctx,
    int num_cells_per_fragment = 10,
    int num_fragments = 1,
    bool overlap = false,
    uint64_t timestamp = 1) {
    // Generate fragments in random order
    std::vector<int> frags(num_fragments);
    std::iota(frags.begin(), frags.end(), 0);
    std::shuffle(frags.begin(), frags.end(), std::random_device{});

    const char* dim_name = "d0";
    const char* attr_name = "a0";

    // Write to SOMAArray
    for (auto i = 0; i < num_fragments; ++i) {
        auto frag_num = frags[i];
        auto soma_array = SOMAArray::open(OpenMode::soma_write, uri, ctx, TimestampRange(timestamp + i, timestamp + i));

        std::vector<int64_t> d0(num_cells_per_fragment);
        for (int j = 0; j < num_cells_per_fragment; j++) {
            // Overlap odd fragments when generating overlaps
            if (overlap && frag_num % 2 == 1) {
                d0[j] = j + num_cells_per_fragment * (frag_num - 1);
            } else {
                d0[j] = j + num_cells_per_fragment * frag_num;
            }
        }
        std::vector<int32_t> a0(num_cells_per_fragment, frag_num);

        // Write data to array
        auto mq = soma_array->create_managed_query("");
        mq.set_layout(ResultOrder::unordered);
        mq.setup_write_column(attr_name, a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.setup_write_column(dim_name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        mq.close();
        soma_array->close();
    }

    // Read from TileDB Array to get expected data
    Array tiledb_array(*ctx->tiledb_ctx(), uri, TILEDB_READ, TemporalPolicy(TimeTravel, timestamp + num_fragments - 1));
    tiledb_array.reopen();

    std::vector<int64_t> expected_d0(num_cells_per_fragment * num_fragments);
    std::vector<int32_t> expected_a0(num_cells_per_fragment * num_fragments);

    Query query(*ctx->tiledb_ctx(), tiledb_array);
    query.set_layout(TILEDB_UNORDERED).set_data_buffer(dim_name, expected_d0).set_data_buffer(attr_name, expected_a0);
    query.submit();

    tiledb_array.close();

    expected_d0.resize(query.result_buffer_elements()[dim_name].second);
    expected_a0.resize(query.result_buffer_elements()[attr_name].second);

    return {expected_d0, expected_a0};
}

};  // namespace

TEST_CASE("SOMAArray: random nnz") {
    auto ctx = std::make_shared<SOMAContext>();

    // Create array
    std::string base_uri = "mem://unit-test-array";
    const char* dim_name = "soma_joinid";
    const char* attr_name = "a0";

    // Create schema
    ArraySchema schema(*ctx->tiledb_ctx(), TILEDB_SPARSE);

    auto dim = Dimension::create<int64_t>(*ctx->tiledb_ctx(), dim_name, {0, std::numeric_limits<int64_t>::max() - 1});

    Domain domain(*ctx->tiledb_ctx());
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int32_t>(*ctx->tiledb_ctx(), attr_name);
    schema.add_attribute(attr);
    schema.set_allows_dups(false);
    schema.check();

    // Create array
    SOMAArray::create(ctx, base_uri, std::move(schema), "NONE", std::nullopt, TimestampRange(0, 2));

    std::vector<size_t> fragment_ids(10);
    std::iota(fragment_ids.begin(), fragment_ids.end(), 0);
    std::shuffle(fragment_ids.begin(), fragment_ids.end(), std::random_device{});

    // Generate 10 random fragments with no overlap
    for (size_t i : fragment_ids) {
        std::vector<int64_t> d0(128 - i);
        std::vector<int32_t> a0(128 - i, i);
        std::iota(d0.begin(), d0.end(), i * 128);

        auto soma_array = SOMAArray::open(OpenMode::soma_write, base_uri, ctx);

        // Write data to array
        auto mq = soma_array->create_managed_query("");
        mq.set_layout(ResultOrder::unordered);
        mq.setup_write_column(attr_name, a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.setup_write_column(dim_name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        mq.close();
        soma_array->close();
    }

    // Generate a randomly overlapping tile
    {
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist(256, 1023);

        auto limit0 = dist(rng);
        auto limit1 = dist(rng);

        auto min = std::min(limit0, limit1);
        auto max = std::max(limit0, limit1);

        std::vector<int64_t> d0(max - min);
        std::vector<int32_t> a0(max - min, max);
        std::iota(d0.begin(), d0.end(), min);

        auto soma_array = SOMAArray::open(OpenMode::soma_write, base_uri, ctx);

        // Write data to array
        auto mq = soma_array->create_managed_query("");
        mq.set_layout(ResultOrder::unordered);
        mq.setup_write_column(attr_name, a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.setup_write_column(dim_name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        mq.close();
        soma_array->close();
    }

    // Get total cell num
    {
        auto soma_array = SOMAArray::open(OpenMode::soma_read, base_uri, ctx);

        uint64_t nnz = soma_array->nnz();
        auto mq = soma_array->create_managed_query();
        mq.set_layout(ResultOrder::unordered);

        uint64_t total_cell_num = 0;
        while (auto batch = mq.read_next()) {
            total_cell_num += (*batch)->num_rows();
        }

        REQUIRE(nnz == total_cell_num);
    }
}

TEST_CASE("SOMAArray: nnz") {
    auto num_fragments = GENERATE(1, 10);
    auto overlap = GENERATE(false, true);
    auto allow_duplicates = true;
    int num_cells_per_fragment = 128;
    auto timestamp = 10;

    const char* dim_name = "d0";
    const char* attr_name = "a0";

    // TODO this use to be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- fragments=" << num_fragments << ", overlap" << overlap << ", allow_duplicates=" << allow_duplicates;

    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();

        // Create array
        std::string base_uri = "mem://unit-test-array";
        auto [uri, expected_nnz] = create_array(
            base_uri, ctx, num_cells_per_fragment, num_fragments, overlap, allow_duplicates);

        // Write at timestamp 10
        auto [expected_d0, expected_a0] = write_array(
            uri, ctx, num_cells_per_fragment, num_fragments, overlap, timestamp);

        // Get total cell num
        auto soma_array = SOMAArray::open(
            OpenMode::soma_read, uri, ctx, TimestampRange(timestamp, timestamp + num_fragments - 1));

        uint64_t nnz = soma_array->nnz();
        REQUIRE(nnz == expected_nnz);
        REQUIRE(nnz >= soma_array->fragment_cell_count());

        std::vector<int64_t> shape = soma_array->shape();
        REQUIRE(shape.size() == 1);
        REQUIRE(shape[0] == std::numeric_limits<int64_t>::max());

        // Check that data from SOMAArray::read_next matches expected data
        auto mq = soma_array->create_managed_query();
        mq.set_layout(ResultOrder::unordered);
        while (!mq.results_complete()) {
            auto arrbuf = mq.read_next().value();
            REQUIRE(arrbuf->names() == std::vector<std::string>({dim_name, attr_name}));
            REQUIRE(arrbuf->num_rows() == nnz);

            auto d0span = arrbuf->at(dim_name)->data<int64_t>();
            auto a0span = arrbuf->at(attr_name)->data<int32_t>();

            std::vector<int64_t> d0col(d0span.begin(), d0span.end());
            std::vector<int32_t> a0col(a0span.begin(), a0span.end());

            REQUIRE(d0col == expected_d0);
            REQUIRE(a0col == expected_a0);
        }
        mq.close();
        soma_array->close();
    }
}

TEST_CASE("SOMAArray: nnz with timestamp") {
    auto num_fragments = GENERATE(1, 10);
    auto overlap = GENERATE(false, true);
    auto allow_duplicates = true;
    int num_cells_per_fragment = 128;

    // TODO this use to be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- fragments=" << num_fragments << ", overlap" << overlap << ", allow_duplicates=" << allow_duplicates;

    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();

        // Create array
        std::string base_uri = "mem://unit-test-array";
        const auto& [uri, expected_nnz] = create_array(
            base_uri, ctx, num_cells_per_fragment, num_fragments, overlap, allow_duplicates);

        // Write at timestamp 10
        write_array(uri, ctx, num_cells_per_fragment, num_fragments, overlap, 10);

        // Write more data to the array at timestamp 40, which will be
        // not be included in the nnz call with a timestamp
        write_array(uri, ctx, num_cells_per_fragment, num_fragments, overlap, 40);

        // Get total cell num at timestamp (0, 20)
        TimestampRange timestamp{0, 20};
        auto soma_array = SOMAArray::open(OpenMode::soma_read, uri, ctx, timestamp);

        uint64_t nnz = soma_array->nnz();
        REQUIRE(nnz == expected_nnz);
        REQUIRE(nnz >= soma_array->fragment_cell_count());
    }
}

TEST_CASE("SOMAArray: nnz with consolidation") {
    auto num_fragments = GENERATE(1, 10);
    auto overlap = GENERATE(false, true);
    auto allow_duplicates = true;
    auto vacuum = GENERATE(false, true);
    int num_cells_per_fragment = 128;

    // TODO this use to be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- fragments=" << num_fragments << ", overlap" << overlap << ", allow_duplicates=" << allow_duplicates;

    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();

        // Create array
        std::string base_uri = "mem://unit-test-array";
        const auto& [uri, expected_nnz] = create_array(
            base_uri, ctx, num_cells_per_fragment, num_fragments, overlap, allow_duplicates);

        // Write at timestamp 10
        write_array(uri, ctx, num_cells_per_fragment, num_fragments, overlap, 10);

        // Write more data to the array at timestamp 20, which will be
        // duplicates of the data written at timestamp 10
        // The duplicates get merged into one fragment during consolidation.
        write_array(uri, ctx, num_cells_per_fragment, num_fragments, overlap, 20);

        // Consolidate and optionally vacuum
        Array::consolidate(*ctx->tiledb_ctx(), uri);
        if (vacuum) {
            Array::vacuum(*ctx->tiledb_ctx(), uri);
        }

        // Get total cell num
        auto soma_array = SOMAArray::open(OpenMode::soma_read, uri, ctx);

        uint64_t nnz = soma_array->nnz();
        if (allow_duplicates) {
            // Since we wrote twice
            REQUIRE(nnz == 2 * expected_nnz);
        } else {
            REQUIRE(nnz == expected_nnz);
        }
    }
}

TEST_CASE("SOMAArray: metadata") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-array";
    const auto& [uri, expected_nnz] = create_array(base_uri, ctx);

    auto soma_array = SOMAArray::open(OpenMode::soma_write, uri, ctx, TimestampRange(1, 1));

    int32_t val = 100;
    soma_array->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_array->close();

    // Read metadata
    soma_array->open(OpenMode::soma_read, TimestampRange(0, 2));
    REQUIRE(soma_array->metadata_num() == 3);
    REQUIRE(soma_array->has_metadata("soma_object_type"));
    REQUIRE(soma_array->has_metadata("soma_encoding_version"));
    REQUIRE(soma_array->has_metadata("md"));
    auto mdval = soma_array->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_array->close();

    // md should not be available at (2, 2)
    soma_array->open(OpenMode::soma_read, TimestampRange(2, 2));
    REQUIRE(soma_array->metadata_num() == 2);
    REQUIRE(soma_array->has_metadata("soma_object_type"));
    REQUIRE(soma_array->has_metadata("soma_encoding_version"));
    REQUIRE(!soma_array->has_metadata("md"));
    soma_array->close();

    // Metadata should also be retrievable in write mode
    soma_array->open(OpenMode::soma_write, TimestampRange(0, 2));
    REQUIRE(soma_array->metadata_num() == 3);
    REQUIRE(soma_array->has_metadata("soma_object_type"));
    REQUIRE(soma_array->has_metadata("soma_encoding_version"));
    REQUIRE(soma_array->has_metadata("md"));
    mdval = soma_array->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    soma_array->delete_metadata("md");
    mdval = soma_array->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_array->close();

    // Confirm delete in read mode
    soma_array->open(OpenMode::soma_read, TimestampRange(0, 2));
    REQUIRE(!soma_array->has_metadata("md"));
    REQUIRE(soma_array->metadata_num() == 2);
}

TEST_CASE("SOMAArray: Test buffer size") {
    // Test soma.init_buffer_bytes by making buffer small
    // enough to read one byte at a time so that read_next
    // must be called 10 times instead of placing all data
    // in buffer within a single read
    std::map<std::string, std::string> cfg;
    cfg["soma.init_buffer_bytes"] = "8";
    auto ctx = std::make_shared<SOMAContext>(cfg);
    REQUIRE(ctx->tiledb_config()["soma.init_buffer_bytes"] == "8");

    std::string base_uri = "mem://unit-test-array";
    auto [uri, expected_nnz] = create_array(base_uri, ctx);
    auto [expected_d0, expected_a0] = write_array(uri, ctx);
    auto soma_array = SOMAArray::open(OpenMode::soma_read, uri, ctx);
    auto mq = soma_array->create_managed_query();

    size_t loops = 0;
    while (auto batch = mq.read_next())
        ++loops;
    REQUIRE(loops == 10);
    soma_array->close();
}

TEST_CASE("SOMAArray: ResultOrder") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-array-result-order";
    auto [uri, expected_nnz] = create_array(base_uri, ctx);
    auto [expected_d0, expected_a0] = write_array(uri, ctx);
    auto soma_array = SOMAArray::open(OpenMode::soma_read, uri, ctx);
    auto mq = soma_array->create_managed_query();
    REQUIRE(mq.result_order() == ResultOrder::automatic);
    mq.set_layout(ResultOrder::rowmajor);
    REQUIRE(mq.result_order() == ResultOrder::rowmajor);
    mq.set_layout(ResultOrder::colmajor);
    REQUIRE(mq.result_order() == ResultOrder::colmajor);
    REQUIRE_THROWS_AS(mq.set_layout(static_cast<ResultOrder>(10)), std::invalid_argument);
}

TEST_CASE("SOMAArray: Write and read back Boolean") {
    std::string uri = "mem://unit-test-array";

    auto ctx = std::make_shared<SOMAContext>();

    const char* dim_name = "d0";
    const char* attr_name = "a0";

    ArraySchema schema(*ctx->tiledb_ctx(), TILEDB_SPARSE);
    auto dim = Dimension::create<int64_t>(*ctx->tiledb_ctx(), dim_name, {0, 7});
    Domain domain(*ctx->tiledb_ctx());
    domain.add_dimension(dim);
    schema.set_domain(domain);
    auto attr = Attribute::create<bool>(*ctx->tiledb_ctx(), attr_name);
    schema.add_attribute(attr);
    schema.set_allows_dups(true);

    SOMAArray::create(ctx, uri, std::move(schema), "NONE");
    auto soma_array = SOMAArray::open(OpenMode::soma_write, uri, ctx);

    auto arrow_schema = make_managed_unique<ArrowSchema>();
    arrow_schema->format = strdup("+s");
    arrow_schema->n_children = 2;
    arrow_schema->dictionary = nullptr;
    arrow_schema->metadata = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->children = (ArrowSchema**)malloc(arrow_schema->n_children * sizeof(ArrowSchema*));
    ArrowSchema* arrow_dim = arrow_schema->children[0] = (ArrowSchema*)malloc(sizeof(ArrowSchema));
    arrow_dim->format = strdup(ArrowAdapter::tdb_to_arrow_type(TILEDB_INT64).c_str());
    arrow_dim->name = strdup(dim_name);
    arrow_dim->n_children = 0;
    arrow_dim->dictionary = nullptr;
    arrow_dim->metadata = nullptr;
    arrow_dim->children = nullptr;
    arrow_dim->release = &ArrowAdapter::release_schema;
    ArrowSchema* arrow_att = arrow_schema->children[1] = (ArrowSchema*)malloc(sizeof(ArrowSchema));
    arrow_att->format = strdup("b");
    arrow_att->name = strdup(attr_name);
    arrow_att->n_children = 0;
    arrow_att->dictionary = nullptr;
    arrow_att->metadata = nullptr;
    arrow_att->children = nullptr;
    arrow_att->release = &ArrowAdapter::release_schema;

    auto arrow_array = make_managed_unique<ArrowArray>();
    arrow_array->length = 0;
    arrow_array->null_count = 0;
    arrow_array->offset = 0;
    arrow_array->n_buffers = 0;
    arrow_array->buffers = nullptr;
    arrow_array->n_children = 2;
    arrow_array->release = &ArrowAdapter::release_array;
    arrow_array->private_data = nullptr;
    arrow_array->children = (ArrowArray**)malloc(2 * sizeof(ArrowArray*));
    arrow_array->dictionary = nullptr;

    auto d0_expected = arrow_array->children[0] = (ArrowArray*)malloc(sizeof(ArrowArray));
    d0_expected->length = 8;
    d0_expected->null_count = 0;
    d0_expected->offset = 0;
    d0_expected->n_buffers = 2;
    d0_expected->release = &ArrowAdapter::release_array;
    d0_expected->private_data = nullptr;
    d0_expected->buffers = (const void**)malloc(2 * sizeof(void*));
    d0_expected->buffers[0] = nullptr;
    d0_expected->buffers[1] = malloc(sizeof(int64_t) * 8);
    d0_expected->n_children = 0;
    d0_expected->children = nullptr;
    d0_expected->dictionary = nullptr;
    int64_t d0_data[] = {0, 1, 2, 3, 4, 5, 6, 7};
    std::memcpy((void*)d0_expected->buffers[1], &d0_data, sizeof(int64_t) * 8);

    auto a0_expected = arrow_array->children[1] = (ArrowArray*)malloc(sizeof(ArrowArray));
    ;
    a0_expected->length = 8;
    a0_expected->null_count = 0;
    a0_expected->offset = 0;
    a0_expected->n_buffers = 2;
    a0_expected->release = &ArrowAdapter::release_array;
    a0_expected->private_data = nullptr;
    a0_expected->buffers = (const void**)malloc(2 * sizeof(void*));
    a0_expected->buffers[0] = nullptr;
    a0_expected->buffers[1] = malloc(sizeof(uint8_t));
    a0_expected->n_children = 0;
    a0_expected->children = nullptr;
    a0_expected->dictionary = nullptr;
    uint8_t a0_data = 0b10101010;
    std::memcpy((void*)a0_expected->buffers[1], &a0_data, sizeof(uint8_t));

    auto mq_write = soma_array->create_managed_query();
    mq_write.set_array_data(arrow_schema.get(), arrow_array.get());
    mq_write.submit_write();
    mq_write.close();
    soma_array->close();

    soma_array = SOMAArray::open(OpenMode::soma_read, uri, ctx);
    auto mq_read = soma_array->create_managed_query();
    auto arrbuf = mq_read.read_next().value();

    auto d0_span = arrbuf->at(dim_name)->data<int64_t>();
    REQUIRE(std::vector<int64_t>(d0_span.begin(), d0_span.end()) == std::vector<int64_t>(d0_data, d0_data + 8));

    auto a0_span = arrbuf->at(attr_name)->data<bool>();
    REQUIRE(
        std::vector<bool>(a0_span.begin(), a0_span.end()) ==
        std::vector<bool>({false, true, false, true, false, true, false, true}));
    soma_array->close();
}
