/**
 * @file   unit_soma_dense_ndarray.cc
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
 * This file manages unit tests for the SOMADenseNDArray class
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
    ArraySchema schema(ctx, TILEDB_DENSE);

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

TEST_CASE("SOMADenseNDArray: basic") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-dense-ndarray-basic";

    SOMADenseNDArray::create(uri, create_schema(*ctx->tiledb_ctx()), ctx);

    auto soma_dense = SOMADenseNDArray::open(uri, OpenMode::read, ctx);
    REQUIRE(soma_dense->uri() == uri);
    REQUIRE(soma_dense->ctx() == ctx);
    REQUIRE(soma_dense->type() == "SOMADenseNDArray");
    REQUIRE(soma_dense->is_sparse() == false);
    auto schema = soma_dense->tiledb_schema();
    REQUIRE(schema->has_attribute("a0"));
    REQUIRE(schema->domain().has_dimension("d0"));
    REQUIRE(soma_dense->ndim() == 1);
    REQUIRE(soma_dense->shape() == std::vector<int64_t>{1001});
    soma_dense->close();

    std::vector<int64_t> d0{1, 10};
    std::vector<int> a0(10, 1);

    soma_dense->open(OpenMode::write);
    soma_dense->set_column_data("a0", a0.size(), a0.data());
    soma_dense->set_column_data("d0", d0.size(), d0.data());
    soma_dense->write();
    soma_dense->close();

    soma_dense->open(OpenMode::read);
    while (auto batch = soma_dense->read_next()) {
        auto arrbuf = batch.value();
        auto d0span = arrbuf->at("d0")->data<int64_t>();
        auto a0span = arrbuf->at("a0")->data<int>();
        REQUIRE(
            std::vector<int64_t>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10} ==
            std::vector<int64_t>(d0span.begin(), d0span.end()));
        REQUIRE(a0 == std::vector<int>(a0span.begin(), a0span.end()));
    }
    soma_dense->close();
}

TEST_CASE("SOMADenseNDArray: metadata") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-dense-ndarray";
    SOMADenseNDArray::create(
        uri, create_schema(*ctx->tiledb_ctx()), ctx, TimestampRange(0, 2));
    auto soma_dense = SOMADenseNDArray::open(
        uri,
        OpenMode::write,
        ctx,
        {},
        ResultOrder::automatic,
        std::pair<uint64_t, uint64_t>(1, 1));

    int32_t val = 100;
    soma_dense->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_dense->close();

    // Read metadata
    soma_dense->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(soma_dense->metadata_num() == 3);
    REQUIRE(soma_dense->has_metadata("soma_object_type"));
    REQUIRE(soma_dense->has_metadata("soma_encoding_version"));
    REQUIRE(soma_dense->has_metadata("md"));
    auto mdval = soma_dense->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_dense->close();

    // md should not be available at (2, 2)
    soma_dense->open(OpenMode::read, TimestampRange(2, 2));
    REQUIRE(soma_dense->metadata_num() == 2);
    REQUIRE(soma_dense->has_metadata("soma_object_type"));
    REQUIRE(soma_dense->has_metadata("soma_encoding_version"));
    REQUIRE(!soma_dense->has_metadata("md"));
    soma_dense->close();

    // Metadata should also be retrievable in write mode
    soma_dense->open(OpenMode::write, TimestampRange(0, 2));
    REQUIRE(soma_dense->metadata_num() == 3);
    REQUIRE(soma_dense->has_metadata("soma_object_type"));
    REQUIRE(soma_dense->has_metadata("soma_encoding_version"));
    REQUIRE(soma_dense->has_metadata("md"));
    mdval = soma_dense->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    soma_dense->delete_metadata("md");
    mdval = soma_dense->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_dense->close();

    // Confirm delete in read mode
    soma_dense->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(!soma_dense->has_metadata("md"));
    REQUIRE(soma_dense->metadata_num() == 2);
}