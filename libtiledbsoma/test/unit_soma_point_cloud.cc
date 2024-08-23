/**
 * @file   unit_soma_point_cloud.cc
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
 * This file manages unit tests for the SOMAPointCloud class
 */

#include "common.h"

#define DIM_MAX 1000

TEST_CASE("SOMAPointCloud: basic") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-point-cloud-basic";

    REQUIRE(!SOMAPointCloud::exists(uri, ctx));

    auto [schema, index_columns] = helper::create_arrow_schema(DIM_MAX);
    SOMAPointCloud::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx);

    REQUIRE(SOMAPointCloud::exists(uri, ctx));
    REQUIRE(!SOMASparseNDArray::exists(uri, ctx));
    REQUIRE(!SOMADenseNDArray::exists(uri, ctx));

    auto soma_point_cloud = SOMAPointCloud::open(uri, OpenMode::read, ctx);
    REQUIRE(soma_point_cloud->uri() == uri);
    REQUIRE(soma_point_cloud->ctx() == ctx);
    REQUIRE(soma_point_cloud->type() == "SOMAPointCloud");
    std::vector<std::string> expected_index_column_names = {"d0"};
    REQUIRE(
        soma_point_cloud->index_column_names() == expected_index_column_names);
    REQUIRE(soma_point_cloud->count() == 0);
    soma_point_cloud->close();
}
