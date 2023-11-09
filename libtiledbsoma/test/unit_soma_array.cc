/**
 * @file   unit_soma_array.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022-2023 TileDB, Inc.
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
    std::shared_ptr<Context> ctx,
    int num_cells_per_fragment = 10,
    int num_fragments = 1,
    bool overlap = false,
    bool allow_duplicates = false) {
    auto vfs = VFS(*ctx);
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    // Create schema
    ArraySchema schema(*ctx, TILEDB_SPARSE);

    auto dim = Dimension::create<int64_t>(
        *ctx, "d0", {0, std::numeric_limits<int64_t>::max() - 1});

    Domain domain(*ctx);
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int>(*ctx, "a0");
    schema.add_attribute(attr);
    schema.set_allows_dups(allow_duplicates);
    schema.check();

    // Create array
    SOMAArray::create(ctx, uri, schema, "NONE");

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

std::tuple<std::vector<int64_t>, std::vector<int>> write_array(
    const std::string& uri,
    std::shared_ptr<Context> ctx,
    int num_cells_per_fragment = 10,
    int num_fragments = 1,
    bool overlap = false,
    uint64_t timestamp = 1) {
    // Generate fragments in random order
    std::vector<int> frags(num_fragments);
    std::iota(frags.begin(), frags.end(), 0);
    std::shuffle(frags.begin(), frags.end(), std::random_device{});

    // Write to SOMAArray
    for (auto i = 0; i < num_fragments; ++i) {
        auto frag_num = frags[i];
        auto soma_array = SOMAArray::open(
            OpenMode::write,
            ctx,
            uri,
            "",
            {},
            "auto",
            ResultOrder::automatic,
            std::pair<uint64_t, uint64_t>(timestamp + i, timestamp + i));

        if (LOG_DEBUG_ENABLED()) {
            soma_array->schema()->dump();
        }

        std::vector<int64_t> d0(num_cells_per_fragment);
        for (int j = 0; j < num_cells_per_fragment; j++) {
            // Overlap odd fragments when generating overlaps
            if (overlap && frag_num % 2 == 1) {
                d0[j] = j + num_cells_per_fragment * (frag_num - 1);
            } else {
                d0[j] = j + num_cells_per_fragment * frag_num;
            }
        }
        std::vector<int> a0(num_cells_per_fragment, frag_num);

        auto array_buffer = std::make_shared<ArrayBuffers>();
        auto tdb_arr = std::make_shared<Array>(*ctx, uri, TILEDB_READ);
        array_buffer->emplace("a0", ColumnBuffer::create(tdb_arr, "a0", a0));
        array_buffer->emplace("d0", ColumnBuffer::create(tdb_arr, "d0", d0));

        // Write data to array
        soma_array->write(array_buffer);
        soma_array->close();
    }

    // Read from TileDB Array to get expected data
    Array tiledb_array(*ctx, uri, TILEDB_READ, timestamp + num_fragments - 1);
    tiledb_array.reopen();

    std::vector<int64_t> expected_d0(num_cells_per_fragment * num_fragments);
    std::vector<int> expected_a0(num_cells_per_fragment * num_fragments);

    Query query(*ctx, tiledb_array);
    query.set_layout(TILEDB_UNORDERED)
        .set_data_buffer("d0", expected_d0)
        .set_data_buffer("a0", expected_a0);
    query.submit();

    tiledb_array.close();

    expected_d0.resize(query.result_buffer_elements()["d0"].second);
    expected_a0.resize(query.result_buffer_elements()["a0"].second);

    return {expected_d0, expected_a0};
}

};  // namespace

TEST_CASE("SOMAArray: nnz") {
    auto num_fragments = GENERATE(1, 10);
    auto overlap = GENERATE(false, true);
    auto allow_duplicates = true;
    int num_cells_per_fragment = 128;
    auto timestamp = 10;

    // TODO this use to be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- fragments=" << num_fragments << ", overlap" << overlap
            << ", allow_duplicates=" << allow_duplicates;

    SECTION(section.str()) {
        auto ctx = std::make_shared<Context>();

        // Create array
        std::string base_uri = "mem://unit-test-array";
        auto [uri, expected_nnz] = create_array(
            base_uri,
            ctx,
            num_cells_per_fragment,
            num_fragments,
            overlap,
            allow_duplicates);

        // Write at timestamp 10
        auto [expected_d0, expected_a0] = write_array(
            uri,
            ctx,
            num_cells_per_fragment,
            num_fragments,
            overlap,
            timestamp);

        // Get total cell num
        auto soma_array = SOMAArray::open(
            OpenMode::read,
            ctx,
            uri,
            "",
            {},
            "auto",
            ResultOrder::automatic,
            std::pair<uint64_t, uint64_t>(
                timestamp, timestamp + num_fragments - 1));

        uint64_t nnz = soma_array->nnz();
        REQUIRE(nnz == expected_nnz);

        std::vector<int64_t> shape = soma_array->shape();
        REQUIRE(shape.size() == 1);
        REQUIRE(shape[0] == std::numeric_limits<int64_t>::max());

        // Check that data from SOMAArray::read_next matches expected data
        while (auto batch = soma_array->read_next()) {
            auto arrbuf = batch.value();
            REQUIRE(arrbuf->names() == std::vector<std::string>({"d0", "a0"}));
            REQUIRE(arrbuf->num_rows() == nnz);

            auto d0span = arrbuf->at("d0")->data<int64_t>();
            auto a0span = arrbuf->at("a0")->data<int>();

            std::vector<int64_t> d0col(d0span.begin(), d0span.end());
            std::vector<int> a0col(a0span.begin(), a0span.end());

            REQUIRE(d0col == expected_d0);
            REQUIRE(a0col == expected_a0);
        }
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
    section << "- fragments=" << num_fragments << ", overlap" << overlap
            << ", allow_duplicates=" << allow_duplicates;

    SECTION(section.str()) {
        auto ctx = std::make_shared<Context>();

        // Create array
        std::string base_uri = "mem://unit-test-array";
        const auto& [uri, expected_nnz] = create_array(
            base_uri,
            ctx,
            num_cells_per_fragment,
            num_fragments,
            overlap,
            allow_duplicates);

        // Write at timestamp 10
        write_array(
            uri, ctx, num_cells_per_fragment, num_fragments, overlap, 10);

        // Write more data to the array at timestamp 40, which will be
        // not be included in the nnz call with a timestamp
        write_array(
            uri, ctx, num_cells_per_fragment, num_fragments, overlap, 40);

        // Get total cell num at timestamp (0, 20)
        std::pair<uint64_t, uint64_t> timestamp{0, 20};
        auto soma_array = SOMAArray::open(
            OpenMode::read,
            ctx,
            uri,
            "nnz",
            {},
            "auto",
            ResultOrder::automatic,
            timestamp);

        uint64_t nnz = soma_array->nnz();
        REQUIRE(nnz == expected_nnz);
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
    section << "- fragments=" << num_fragments << ", overlap" << overlap
            << ", allow_duplicates=" << allow_duplicates;

    SECTION(section.str()) {
        auto ctx = std::make_shared<Context>();

        // Create array
        std::string base_uri = "mem://unit-test-array";
        const auto& [uri, expected_nnz] = create_array(
            base_uri,
            ctx,
            num_cells_per_fragment,
            num_fragments,
            overlap,
            allow_duplicates);

        // Write at timestamp 10
        write_array(
            uri, ctx, num_cells_per_fragment, num_fragments, overlap, 10);

        // Write more data to the array at timestamp 20, which will be
        // duplicates of the data written at timestamp 10
        // The duplicates get merged into one fragment during consolidation.
        write_array(
            uri, ctx, num_cells_per_fragment, num_fragments, overlap, 20);

        // Consolidate and optionally vacuum
        Array::consolidate(*ctx, uri);
        if (vacuum) {
            Array::vacuum(*ctx, uri);
        }

        // Get total cell num
        auto soma_array = SOMAArray::open(
            OpenMode::read,
            ctx,
            uri,
            "nnz",
            {},
            "auto",
            ResultOrder::automatic);

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
    auto ctx = std::make_shared<Context>();

    std::string base_uri = "mem://unit-test-array";
    const auto& [uri, expected_nnz] = create_array(base_uri, ctx);

    auto soma_array = SOMAArray::open(
        OpenMode::write,
        ctx,
        uri,
        "metadata_test",
        {},
        "auto",
        ResultOrder::automatic,
        std::pair<uint64_t, uint64_t>(1, 1));
    int32_t val = 100;
    soma_array->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_array->close();

    soma_array->open(OpenMode::read, std::pair<uint64_t, uint64_t>(1, 1));
    REQUIRE(soma_array->metadata_num() == 2);
    REQUIRE(soma_array->has_metadata("soma_object_type") == true);
    REQUIRE(soma_array->has_metadata("md") == true);

    auto mdval = soma_array->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_array->close();

    soma_array->open(OpenMode::write, std::pair<uint64_t, uint64_t>(2, 2));
    // Metadata should also be retrievable in write mode
    mdval = soma_array->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_array->delete_metadata("md");
    mdval = soma_array->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_array->close();

    soma_array->open(OpenMode::read, std::pair<uint64_t, uint64_t>(3, 3));
    REQUIRE(soma_array->has_metadata("md") == false);
    REQUIRE(soma_array->metadata_num() == 1);
    soma_array->close();
}

TEST_CASE("SOMAArray: Test buffer size") {
    // Test soma.init_buffer_bytes by making buffer small
    // enough to read one byte at a time so that read_next
    // must be called 10 times instead of placing all data
    // in buffer within a single read
    Config cfg;
    cfg["soma.init_buffer_bytes"] = 8;
    auto ctx = std::make_shared<Context>(cfg);

    std::string base_uri = "mem://unit-test-array";
    auto [uri, expected_nnz] = create_array(base_uri, ctx);
    auto [expected_d0, expected_a0] = write_array(uri, ctx);
    auto soma_array = SOMAArray::open(OpenMode::read, ctx, uri);

    size_t loops = 0;
    while (auto batch = soma_array->read_next())
        ++loops;
    REQUIRE(loops == 10);
    soma_array->close();
}

TEST_CASE("SOMAArray: Enumeration") {
    std::string uri = "mem://unit-test-array-enmr";
    auto ctx = std::make_shared<Context>();
    ArraySchema schema(*ctx, TILEDB_SPARSE);

    auto dim = Dimension::create<int64_t>(
        *ctx, "d", {0, std::numeric_limits<int64_t>::max() - 1});

    Domain dom(*ctx);
    dom.add_dimension(dim);
    schema.set_domain(dom);

    std::vector<std::string> vals = {"red", "blue", "green"};
    auto enmr = Enumeration::create(*ctx, "rbg", vals);
    ArraySchemaExperimental::add_enumeration(*ctx, schema, enmr);

    auto attr = Attribute::create<int>(*ctx, "a");
    AttributeExperimental::set_enumeration_name(*ctx, attr, "rbg");
    schema.add_attribute(attr);

    Array::create(uri, schema);

    auto soma_array = SOMAArray::open(OpenMode::read, ctx, uri);
    auto attr_to_enum = soma_array->get_attr_to_enum_mapping();
    REQUIRE(attr_to_enum.size() == 1);
    REQUIRE(attr_to_enum.at("a").name() == "rbg");
    REQUIRE(soma_array->get_enum_label_on_attr("a"));
    REQUIRE(soma_array->attr_has_enum("a"));
}

TEST_CASE("SOMAArray: ResultOrder") {
    auto ctx = std::make_shared<Context>();
    std::string base_uri = "mem://unit-test-array-result-order";
    auto [uri, expected_nnz] = create_array(base_uri, ctx);
    auto [expected_d0, expected_a0] = write_array(uri, ctx);
    auto soma_array = SOMAArray::open(OpenMode::read, ctx, uri);
    REQUIRE(soma_array->result_order() == ResultOrder::automatic);
    soma_array->reset({}, "auto", ResultOrder::rowmajor);
    REQUIRE(soma_array->result_order() == ResultOrder::rowmajor);
    soma_array->reset({}, "auto", ResultOrder::colmajor);
    REQUIRE(soma_array->result_order() == ResultOrder::colmajor);
    REQUIRE_THROWS_AS(
        soma_array->reset({}, "auto", static_cast<ResultOrder>(3)),
        std::invalid_argument);
}
