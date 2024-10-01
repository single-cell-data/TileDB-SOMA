/**
 * @file   unit_soma_point_cloud_dataframe.cc
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
 * This file manages unit tests for the SOMAPointCloudDataFrame class
 */

#include "common.h"

const int64_t SOMA_JOINID_DIM_MAX = 99;

TEST_CASE("SOMAPointCloudDataFrame: basic", "[SOMAPointCloudDataFrame]") {
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();
        std::string uri{"mem://unit-test-point-cloud-basic"};
        PlatformConfig platform_config{};

        std::vector<helper::DimInfo> dim_infos({
            helper::DimInfo(
                {.name = "soma_joinid",
                 .tiledb_datatype = TILEDB_INT64,
                 .dim_max = SOMA_JOINID_DIM_MAX,
                 .string_lo = "N/A",
                 .string_hi = "N/A",
                 .use_current_domain = use_current_domain}),
            helper::DimInfo(
                {.name = "x",
                 .tiledb_datatype = TILEDB_UINT32,
                 .dim_max = 100,
                 .string_lo = "N/A",
                 .string_hi = "N/A",
                 .use_current_domain = use_current_domain}),
            helper::DimInfo(
                {.name = "y",
                 .tiledb_datatype = TILEDB_UINT32,
                 .dim_max = 100,
                 .string_lo = "N/A",
                 .string_hi = "N/A",
                 .use_current_domain = use_current_domain}),
        });

        std::vector<helper::AttrInfo> attr_infos({helper::AttrInfo(
            {.name = "radius", .tiledb_datatype = TILEDB_FLOAT64})});

        // Check the point cloud doesn't exist yet.
        REQUIRE(!SOMAPointCloudDataFrame::exists(uri, ctx));

        // Create the point cloud.
        auto [schema, index_columns] =
            helper::create_arrow_schema_and_index_columns(
                dim_infos, attr_infos);
        SOMAPointCloudDataFrame::create(
            uri,
            std::move(schema),
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            platform_config,
            std::nullopt);

        // Check the point cloud exists and it cannot be read as a different
        // object.
        REQUIRE(SOMAPointCloudDataFrame::exists(uri, ctx));
        REQUIRE(!SOMASparseNDArray::exists(uri, ctx));
        REQUIRE(!SOMADenseNDArray::exists(uri, ctx));
        REQUIRE(!SOMADataFrame::exists(uri, ctx));

        auto soma_point_cloud = SOMAPointCloudDataFrame::open(
            uri,
            OpenMode::read,
            ctx,
            {},  // column_names,
            ResultOrder::automatic,
            std::nullopt);
        REQUIRE(soma_point_cloud->uri() == uri);
        REQUIRE(soma_point_cloud->ctx() == ctx);
        REQUIRE(soma_point_cloud->type() == "SOMAPointCloudDataFrame");
        std::vector<std::string> expected_index_column_names = {
            dim_infos[0].name, dim_infos[1].name, dim_infos[2].name};
        REQUIRE(
            soma_point_cloud->index_column_names() ==
            expected_index_column_names);
        REQUIRE(soma_point_cloud->nnz() == 0);
        soma_point_cloud->close();

        // Create vectors of data for writing.
        std::vector<int64_t> d0(10);
        std::iota(d0.begin(), d0.end(), 0);
        std::vector<uint32_t> d1(10);
        std::iota(d1.begin(), d1.end(), 1);
        std::vector<uint32_t> d2(10, 10);
        std::iota(d2.begin(), d2.end(), 0.0);
        std::vector<double> a0(10, 1.0);

        // Write to point cloud.
        soma_point_cloud = SOMAPointCloudDataFrame::open(
            uri,
            OpenMode::write,
            ctx,
            {},  // column_names
            ResultOrder::automatic,
            std::nullopt);
        soma_point_cloud->set_column_data(
            dim_infos[0].name, d0.size(), d0.data());
        soma_point_cloud->set_column_data(
            dim_infos[1].name, d1.size(), d1.data());
        soma_point_cloud->set_column_data(
            dim_infos[2].name, d2.size(), d2.data());
        soma_point_cloud->set_column_data(
            attr_infos[0].name, a0.size(), a0.data());
        soma_point_cloud->write();
        soma_point_cloud->close();

        // Read back the data.
        soma_point_cloud = SOMAPointCloudDataFrame::open(
            uri,
            OpenMode::read,
            ctx,
            {},  // column_names,
            ResultOrder::automatic,
            std::nullopt);
        while (auto batch = soma_point_cloud->read_next()) {
            auto arrbuf = batch.value();
            auto d0span = arrbuf->at(dim_infos[0].name)->data<int64_t>();
            auto d1span = arrbuf->at(dim_infos[1].name)->data<uint32_t>();
            auto d2span = arrbuf->at(dim_infos[2].name)->data<uint32_t>();
            auto a0span = arrbuf->at(attr_infos[0].name)->data<double>();
            CHECK(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
            CHECK(d1 == std::vector<uint32_t>(d1span.begin(), d1span.end()));
            CHECK(d2 == std::vector<uint32_t>(d2span.begin(), d2span.end()));
            CHECK(a0 == std::vector<double>(a0span.begin(), a0span.end()));
        }
        soma_point_cloud->close();

        auto soma_object = SOMAObject::open(uri, OpenMode::read, ctx);
        REQUIRE(soma_object->uri() == uri);
        REQUIRE(soma_object->type() == "SOMAPointCloudDataFrame");
        soma_object->close();
    }
}
