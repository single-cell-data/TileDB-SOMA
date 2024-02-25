/**
 * @file   unit_soma_dataframe.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
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
 * This file manages unit tests for the SOMADataFrame class
 */

#include "common.h"

TEST_CASE("SOMADataFrame: basic") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-dataframe-basic";

    auto [schema, index_columns] = helper::create_arrow_schema();
    SOMADataFrame::create(uri, schema, index_columns, ctx);

    auto soma_dataframe = SOMADataFrame::open(uri, OpenMode::read, ctx);
    REQUIRE(soma_dataframe->uri() == uri);
    REQUIRE(soma_dataframe->ctx() == ctx);
    REQUIRE(soma_dataframe->type() == "SOMADataFrame");
    std::vector<std::string> expected_index_column_names = {"d0"};
    REQUIRE(
        soma_dataframe->index_column_names() == expected_index_column_names);
    REQUIRE(soma_dataframe->count() == 0);
    soma_dataframe->close();

    std::vector<int64_t> d0(10);
    for (int j = 0; j < 10; j++)
        d0[j] = j;
    std::vector<int> a0(10, 1);

    auto array_buffer = std::make_shared<ArrayBuffers>();
    auto tdb_arr = std::make_shared<Array>(
        *ctx->tiledb_ctx(), uri, TILEDB_READ);
    auto col_a0 = ColumnBuffer::create(tdb_arr, "a0");
    auto col_d0 = ColumnBuffer::create(tdb_arr, "d0");
    col_a0->set_data(a0);
    col_d0->set_data(d0);
    array_buffer->emplace("a0", col_a0);
    array_buffer->emplace("d0", col_d0);

    soma_dataframe = SOMADataFrame::open(uri, OpenMode::write, ctx);
    soma_dataframe->write(array_buffer);
    soma_dataframe->close();

    soma_dataframe = SOMADataFrame::open(uri, OpenMode::read, ctx);
    while (auto batch = soma_dataframe->read_next()) {
        auto arrbuf = batch.value();
        auto d0span = arrbuf->at("d0")->data<int64_t>();
        auto a0span = arrbuf->at("a0")->data<int>();
        REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
        REQUIRE(a0 == std::vector<int>(a0span.begin(), a0span.end()));
    }
    soma_dataframe->close();

    auto soma_object = SOMAObject::open(uri, OpenMode::read, ctx);
    REQUIRE(soma_object->uri() == uri);
    REQUIRE(soma_object->type() == "SOMADataFrame");
    soma_object->close();
}

TEST_CASE("SOMADataFrame: metadata") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-collection";
    auto [schema, index_columns] = helper::create_arrow_schema();
    SOMADataFrame::create(uri, schema, index_columns, ctx);
    auto soma_dataframe = SOMADataFrame::open(
        uri,
        OpenMode::write,
        ctx,
        {},
        ResultOrder::automatic,
        std::pair<uint64_t, uint64_t>(1, 1));
    int32_t val = 100;
    soma_dataframe->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_dataframe->close();

    soma_dataframe->open(OpenMode::read, std::pair<uint64_t, uint64_t>(1, 1));
    REQUIRE(soma_dataframe->metadata_num() == 2);
    REQUIRE(soma_dataframe->has_metadata("soma_object_type") == true);
    REQUIRE(soma_dataframe->has_metadata("md") == true);

    auto mdval = soma_dataframe->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_dataframe->close();

    soma_dataframe->open(OpenMode::write, std::pair<uint64_t, uint64_t>(2, 2));
    // Metadata should also be retrievable in write mode
    mdval = soma_dataframe->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_dataframe->delete_metadata("md");
    mdval = soma_dataframe->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_dataframe->close();

    soma_dataframe->open(OpenMode::read, std::pair<uint64_t, uint64_t>(3, 3));
    REQUIRE(soma_dataframe->has_metadata("md") == false);
    REQUIRE(soma_dataframe->metadata_num() == 1);
    soma_dataframe->close();
}