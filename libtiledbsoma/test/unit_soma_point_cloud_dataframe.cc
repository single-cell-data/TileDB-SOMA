/**
 * @file   unit_soma_point_cloud_dataframe.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMAPointCloudDataFrame class
 */

#include <format>
#include "common.h"

const int64_t SOMA_JOINID_DIM_MAX = 99;

TEST_CASE(
    "SOMAPointCloudDataFrame: basic", "[point_cloud_dataframe][spatial]") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri{"mem://unit-test-point-cloud-basic"};
    PlatformConfig platform_config{};

    std::vector<helper::DimInfo> dim_infos({
        helper::DimInfo(
            {.name = "soma_joinid",
             .tiledb_datatype = TILEDB_INT64,
             .dim_max = SOMA_JOINID_DIM_MAX,
             .string_lo = "N/A",
             .string_hi = "N/A"}),
        helper::DimInfo(
            {.name = "x",
             .tiledb_datatype = TILEDB_UINT32,
             .dim_max = 100,
             .string_lo = "N/A",
             .string_hi = "N/A"}),
        helper::DimInfo(
            {.name = "y",
             .tiledb_datatype = TILEDB_UINT32,
             .dim_max = 100,
             .string_lo = "N/A",
             .string_hi = "N/A"}),
    });

    std::vector<helper::AttrInfo> attr_infos({helper::AttrInfo(
        {.name = "radius", .tiledb_datatype = TILEDB_FLOAT64})});

    // Check the point cloud doesn't exist yet.
    REQUIRE(!SOMAPointCloudDataFrame::exists(uri, ctx));

    // Create the point cloud.
    auto [schema, index_columns] =
        helper::create_arrow_schema_and_index_columns(dim_infos, attr_infos);
    SOMACoordinateSpace coord_space{};
    SOMAPointCloudDataFrame::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        coord_space,
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
        uri, OpenMode::read, ctx, std::nullopt);
    REQUIRE(soma_point_cloud->uri() == uri);
    REQUIRE(soma_point_cloud->ctx() == ctx);
    REQUIRE(soma_point_cloud->type() == "SOMAPointCloudDataFrame");
    std::vector<std::string> expected_index_column_names = {
        dim_infos[0].name, dim_infos[1].name, dim_infos[2].name};
    REQUIRE(
        soma_point_cloud->index_column_names() == expected_index_column_names);
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
    soma_point_cloud = SOMAPointCloudDataFrame::open(uri, OpenMode::write, ctx);
    auto mq = ManagedQuery(*soma_point_cloud, ctx->tiledb_ctx());
    mq.setup_write_column(
        dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
    mq.setup_write_column(
        dim_infos[1].name, d1.size(), d1.data(), (uint64_t*)nullptr);
    mq.setup_write_column(
        dim_infos[2].name, d2.size(), d2.data(), (uint64_t*)nullptr);
    mq.setup_write_column(
        attr_infos[0].name, a0.size(), a0.data(), (uint64_t*)nullptr);
    mq.submit_write();
    soma_point_cloud->close();

    // Read back the data.
    soma_point_cloud = SOMAPointCloudDataFrame::open(
        uri, OpenMode::read, ctx, std::nullopt);
    mq = ManagedQuery(*soma_point_cloud, ctx->tiledb_ctx());
    while (auto batch = mq.read_next()) {
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
    CHECK(soma_point_cloud->has_metadata("soma_encoding_version"));
    CHECK(soma_point_cloud->has_metadata("soma_spatial_encoding_version"));
    auto point_cloud_coord_space = soma_point_cloud->coordinate_space();
    CHECK(point_cloud_coord_space == coord_space);
    soma_point_cloud->close();

    auto soma_object = SOMAObject::open(uri, OpenMode::read, ctx);
    REQUIRE(soma_object->uri() == uri);
    REQUIRE(soma_object->type() == "SOMAPointCloudDataFrame");
    soma_object->close();
}
