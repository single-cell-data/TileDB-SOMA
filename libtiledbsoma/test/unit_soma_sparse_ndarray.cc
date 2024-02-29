/**
 * @file   unit_soma_sparse_ndarray.cc
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
 * This file manages unit tests for the SOMASparseNDArray class
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
ArraySchema create_schema(Context& ctx, bool allow_duplicates = false) {
    // Create schema
    ArraySchema schema(ctx, TILEDB_SPARSE);

    auto dim = Dimension::create<int64_t>(ctx, "d0", {0, 1000});

    Domain domain(ctx);
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int>(ctx, "a0");
    schema.add_attribute(attr);
    schema.set_allows_dups(allow_duplicates);
    schema.check();

    return schema;
}
};  // namespace

TEST_CASE("SOMASparseNDArray: basic") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-sparse-ndarray-basic";

    auto soma_sparse = SOMASparseNDArray::create(
        uri, create_schema(*ctx->tiledb_ctx()), ctx);

    std::vector<int64_t> d0(10);
    for (int j = 0; j < 10; j++)
        d0[j] = j;
    std::vector<int> a0(10, 1);

    auto array_buffer = std::make_shared<ArrayBuffers>();
    auto tdb_arr = std::make_shared<Array>(
        *ctx->tiledb_ctx(), uri, TILEDB_READ);
    array_buffer->emplace("a0", ColumnBuffer::create(tdb_arr, "a0", a0));
    array_buffer->emplace("d0", ColumnBuffer::create(tdb_arr, "d0", d0));

    soma_sparse->write(array_buffer);
    soma_sparse->close();

    soma_sparse->open(OpenMode::read);
    REQUIRE(soma_sparse->uri() == uri);
    REQUIRE(soma_sparse->ctx() == ctx);
    REQUIRE(soma_sparse->type() == "SOMASparseNDArray");
    while (auto batch = soma_sparse->read_next()) {
        auto arrbuf = batch.value();
        auto d0span = arrbuf->at("d0")->data<int64_t>();
        auto a0span = arrbuf->at("a0")->data<int>();
        REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
        REQUIRE(a0 == std::vector<int>(a0span.begin(), a0span.end()));
    }
    soma_sparse->close();
}

TEST_CASE("SOMASparseNDArray: metadata") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-sparse-ndarray";
    SOMASparseNDArray::create(
        uri, create_schema(*ctx->tiledb_ctx()), ctx, TimestampRange(0, 2));
    auto soma_sparse = SOMASparseNDArray::open(
        uri,
        OpenMode::write,
        ctx,
        {},
        ResultOrder::automatic,
        std::pair<uint64_t, uint64_t>(1, 1));

    int32_t val = 100;
    soma_sparse->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_sparse->close();

    // Read metadata
    soma_sparse->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(soma_sparse->metadata_num() == 3);
    REQUIRE(soma_sparse->has_metadata("soma_object_type"));
    REQUIRE(soma_sparse->has_metadata("soma_encoding_version"));
    REQUIRE(soma_sparse->has_metadata("md"));
    auto mdval = soma_sparse->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_sparse->close();

    // md should not be available at (2, 2)
    soma_sparse->open(OpenMode::read, TimestampRange(2, 2));
    REQUIRE(soma_sparse->metadata_num() == 2);
    REQUIRE(soma_sparse->has_metadata("soma_object_type"));
    REQUIRE(soma_sparse->has_metadata("soma_encoding_version"));
    REQUIRE(!soma_sparse->has_metadata("md"));
    soma_sparse->close();

    // Metadata should also be retrievable in write mode
    soma_sparse->open(OpenMode::write, TimestampRange(0, 2));
    REQUIRE(soma_sparse->metadata_num() == 3);
    REQUIRE(soma_sparse->has_metadata("soma_object_type"));
    REQUIRE(soma_sparse->has_metadata("soma_encoding_version"));
    REQUIRE(soma_sparse->has_metadata("md"));
    mdval = soma_sparse->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    soma_sparse->delete_metadata("md");
    mdval = soma_sparse->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_sparse->close();

    // Confirm delete in read mode
    soma_sparse->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(!soma_sparse->has_metadata("md"));
    REQUIRE(soma_sparse->metadata_num() == 2);
}