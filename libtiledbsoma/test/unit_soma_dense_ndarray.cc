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

#include "common.h"

TEST_CASE("SOMADenseNDArray: basic") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-dense-ndarray-basic";

    auto index_columns = helper::create_column_index_info();
    SOMADenseNDArray::create(
        uri,
        "l",
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        PlatformConfig(),
        TimestampRange(0, 2));

    auto soma_dense = SOMADenseNDArray::open(uri, OpenMode::read, ctx);
    REQUIRE(soma_dense->uri() == uri);
    REQUIRE(soma_dense->ctx() == ctx);
    REQUIRE(soma_dense->type() == "SOMADenseNDArray");
    REQUIRE(soma_dense->is_sparse() == false);
    auto schema = soma_dense->tiledb_schema();
    REQUIRE(schema->has_attribute("soma_data"));
    REQUIRE(schema->array_type() == TILEDB_DENSE);
    REQUIRE(schema->domain().has_dimension("soma_dim_0"));
    REQUIRE(soma_dense->ndim() == 1);
    REQUIRE(soma_dense->shape() == std::vector<int64_t>{1001});
    soma_dense->close();

    std::vector<int64_t> d0{1, 10};
    std::vector<int> a0(10, 1);

    soma_dense->open(OpenMode::write);
    soma_dense->set_column_data("soma_data", a0.size(), a0.data());
    soma_dense->set_column_data("soma_dim_0", d0.size(), d0.data());
    soma_dense->write();
    soma_dense->close();

    soma_dense->open(OpenMode::read);
    while (auto batch = soma_dense->read_next()) {
        auto arrbuf = batch.value();
        auto a0span = arrbuf->at("soma_data")->data<int>();
        REQUIRE(a0 == std::vector<int>(a0span.begin(), a0span.end()));
    }
    soma_dense->close();
}

TEST_CASE("SOMADenseNDArray: metadata") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-dense-ndarray";

    auto index_columns = helper::create_column_index_info();
    SOMASparseNDArray::create(
        uri,
        "l",
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        PlatformConfig(),
        TimestampRange(0, 2));

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
