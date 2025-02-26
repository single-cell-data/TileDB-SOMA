/**
 * @file   unit_soma_dense_ndarray.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMADenseNDArray class
 */

#include <format>
#include "common.h"

TEST_CASE("SOMADenseNDArray: basic", "[SOMADenseNDArray]") {
    // Core uses domain & current domain like (0, 999); SOMA uses shape like
    // 1000. We want to carefully and explicitly test here that there aren't any
    // off-by-one errors.
    int64_t dim_max = 999;

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-dense-ndarray-basic";
    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        attr_tiledb_datatype);

    REQUIRE(!SOMADenseNDArray::exists(uri, ctx));

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    if (helper::have_dense_current_domain_support()) {
        SOMADenseNDArray::create(
            uri,
            dim_arrow_format,
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            PlatformConfig(),
            TimestampRange(0, 2));

        auto dnda = SOMADenseNDArray::open(uri, OpenMode::read, ctx);
        REQUIRE(dnda->shape() == std::vector<int64_t>{dim_max + 1});
        dnda->close();
    } else {
        REQUIRE_THROWS(SOMADenseNDArray::create(
            uri,
            dim_arrow_format,
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            PlatformConfig(),
            TimestampRange(0, 2)));
    }
}

TEST_CASE("SOMADenseNDArray: platform_config", "[SOMADenseNDArray]") {
    int64_t dim_max = 999;
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-dense-ndarray-platform-config";
    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t tiledb_datatype = TILEDB_INT64;
    std::string arrow_format = ArrowAdapter::tdb_to_arrow_type(tiledb_datatype);

    PlatformConfig platform_config;
    platform_config.dense_nd_array_dim_zstd_level = 6;

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    if (helper::have_dense_current_domain_support()) {
        SOMADenseNDArray::create(
            uri,
            arrow_format,
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            platform_config);

        auto dnda = SOMADenseNDArray::open(uri, OpenMode::read, ctx);
        auto dim_filter = dnda->tiledb_schema()
                              ->domain()
                              .dimension(dim_name)
                              .filter_list()
                              .filter(0);
        REQUIRE(dim_filter.filter_type() == TILEDB_FILTER_ZSTD);
        REQUIRE(dim_filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL) == 6);

        dnda->close();

    } else {
        REQUIRE_THROWS(SOMADenseNDArray::create(
            uri,
            arrow_format,
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            platform_config));
    }
}

TEST_CASE("SOMADenseNDArray: metadata", "[SOMADenseNDArray]") {
    int64_t dim_max = 999;
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-dense-ndarray";
    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t tiledb_datatype = TILEDB_INT64;
    std::string arrow_format = ArrowAdapter::tdb_to_arrow_type(tiledb_datatype);

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMADenseNDArray::create(
        uri,
        arrow_format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        PlatformConfig(),
        TimestampRange(0, 1));

    // TO DO: do more data writes and readbacks here in C++ tests.
    // https://github.com/single-cell-data/TileDB-SOMA/issues/3721
    auto dnda = SOMADenseNDArray::open(
        uri, OpenMode::write, ctx, TimestampRange(0, 2));

    int32_t val = 100;
    dnda->set_metadata("md", TILEDB_INT32, 1, &val);
    dnda->close();

    // Read metadata
    dnda->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(dnda->metadata_num() == 3);
    REQUIRE(dnda->has_metadata("soma_object_type"));
    REQUIRE(dnda->has_metadata("soma_encoding_version"));
    REQUIRE(dnda->has_metadata("md"));
    auto mdval = dnda->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    dnda->close();

    // md should not be available at (0, 1)
    dnda->open(OpenMode::read, TimestampRange(0, 1));
    REQUIRE(dnda->metadata_num() == 2);
    REQUIRE(dnda->has_metadata("soma_object_type"));
    REQUIRE(dnda->has_metadata("soma_encoding_version"));
    REQUIRE(!dnda->has_metadata("md"));
    dnda->close();

    // Metadata should also be retrievable in write mode
    dnda->open(OpenMode::write);
    REQUIRE(dnda->metadata_num() == 3);
    REQUIRE(dnda->has_metadata("soma_object_type"));
    REQUIRE(dnda->has_metadata("soma_encoding_version"));
    REQUIRE(dnda->has_metadata("md"));
    mdval = dnda->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    dnda->delete_metadata("md");
    mdval = dnda->get_metadata("md");
    REQUIRE(!mdval.has_value());
    dnda->close();

    // Confirm delete in read mode
    dnda->open(OpenMode::read);
    REQUIRE(!dnda->has_metadata("md"));
    REQUIRE(dnda->metadata_num() == 2);
}
