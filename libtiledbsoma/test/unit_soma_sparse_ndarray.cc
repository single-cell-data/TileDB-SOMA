/**
 * @file   unit_soma_sparse_ndarray.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMASparseNDArray class
 */

#include <format>
#include "common.h"

TEST_CASE("SOMASparseNDArray: basic", "[SOMASparseNDArray]") {
    // Core uses domain & current domain like (0, 999); SOMA uses shape like
    // 1000. We want to carefully and explicitly test here that there aren't any
    // off-by-one errors.
    int64_t dim_max = 999;
    int64_t shape = 1000;

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-sparse-ndarray-basic";
    std::string dim_name = "soma_dim_0";
    std::string attr_name = "soma_data";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        attr_tiledb_datatype);

    REQUIRE(!SOMASparseNDArray::exists(uri, ctx));

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(
        uri,
        attr_arrow_format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        PlatformConfig(),
        TimestampRange(0, 2));

    REQUIRE(SOMASparseNDArray::exists(uri, ctx));
    REQUIRE(!SOMADataFrame::exists(uri, ctx));
    REQUIRE(!SOMADenseNDArray::exists(uri, ctx));

    auto snda = SOMASparseNDArray::open(uri, OpenMode::read, ctx);
    REQUIRE(snda->uri() == uri);
    REQUIRE(snda->ctx() == ctx);
    REQUIRE(snda->type() == "SOMASparseNDArray");
    REQUIRE(snda->is_sparse() == true);
    REQUIRE(snda->soma_data_type() == attr_arrow_format);
    auto schema = snda->tiledb_schema();
    REQUIRE(schema->has_attribute(attr_name));
    REQUIRE(schema->array_type() == TILEDB_SPARSE);
    REQUIRE(schema->domain().has_dimension(dim_name));
    REQUIRE(snda->ndim() == 1);
    REQUIRE(snda->nnz() == 0);

    auto expect = std::vector<int64_t>({shape});
    REQUIRE(snda->shape() == expect);

    snda->close();

    std::vector<int64_t> d0(10);
    for (int j = 0; j < 10; j++)
        d0[j] = j;
    std::vector<int32_t> a0(10, 1);

    // A write in read mode should fail
    snda->open(OpenMode::read);
    auto mq = ManagedQuery(*snda, ctx->tiledb_ctx());
    mq.setup_write_column(dim_name, d0.size(), d0.data(), (uint64_t*)nullptr);
    mq.setup_write_column(attr_name, a0.size(), a0.data(), (uint64_t*)nullptr);
    REQUIRE_THROWS(mq.submit_write());
    snda->close();

    snda->open(OpenMode::write);
    mq = ManagedQuery(*snda, ctx->tiledb_ctx());
    mq.setup_write_column(dim_name, d0.size(), d0.data(), (uint64_t*)nullptr);
    mq.setup_write_column(attr_name, a0.size(), a0.data(), (uint64_t*)nullptr);
    mq.submit_write();
    snda->close();

    snda->open(OpenMode::read);
    mq = ManagedQuery(*snda, ctx->tiledb_ctx());
    while (auto batch = mq.read_next()) {
        auto arrbuf = batch.value();
        auto d0span = arrbuf->at(dim_name)->data<int64_t>();
        auto a0span = arrbuf->at(attr_name)->data<int32_t>();
        REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
        REQUIRE(a0 == std::vector<int32_t>(a0span.begin(), a0span.end()));
    }
    snda->close();

    std::vector<int64_t> d0b({dim_max, dim_max + 1});
    std::vector<int64_t> a0b({30, 40});

    // Try out-of-bounds write before resize.
    // * Without current domain support: this should throw since it's
    //   outside the (immutable) doqain.
    // * With current domain support: this should throw since it's outside
    // the (mutable) current domain.
    snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
    mq = ManagedQuery(*snda, ctx->tiledb_ctx());
    mq.setup_write_column(dim_name, d0b.size(), d0b.data(), (uint64_t*)nullptr);
    mq.setup_write_column(
        attr_name, a0b.size(), a0b.data(), (uint64_t*)nullptr);
    REQUIRE_THROWS(mq.submit_write());
    snda->close();

    auto new_shape = std::vector<int64_t>({shape * 2});

    snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
    // Should throw since this already has a shape (core current
    // domain).
    REQUIRE_THROWS(snda->upgrade_shape(new_shape, "testing"));
    snda->resize(new_shape, "testing");
    snda->close();

    // Try out-of-bounds write after resize.
    snda->open(OpenMode::write);
    mq = ManagedQuery(*snda, ctx->tiledb_ctx());
    mq.setup_write_column(dim_name, d0b.size(), d0b.data(), (uint64_t*)nullptr);
    mq.setup_write_column(
        attr_name, a0b.size(), a0b.data(), (uint64_t*)nullptr);
    // Implicitly checking for no throw
    mq.submit_write();
    snda->close();

    snda->open(OpenMode::read);
    REQUIRE(snda->shape() == new_shape);
    snda->close();
}

TEST_CASE("SOMASparseNDArray: platform_config", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-dataframe-platform-config";
    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        attr_tiledb_datatype);

    PlatformConfig platform_config;
    platform_config.sparse_nd_array_dim_zstd_level = 6;

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(
        uri,
        attr_arrow_format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        platform_config);

    auto soma_dataframe = SOMASparseNDArray::open(uri, OpenMode::read, ctx);
    auto dim_filter = soma_dataframe->tiledb_schema()
                          ->domain()
                          .dimension(dim_name)
                          .filter_list()
                          .filter(0);
    REQUIRE(dim_filter.filter_type() == TILEDB_FILTER_ZSTD);
    REQUIRE(dim_filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL) == 6);

    soma_dataframe->close();
}

TEST_CASE("SOMASparseNDArray: metadata", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-sparse-ndarray";
    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        attr_tiledb_datatype);

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(
        uri,
        attr_arrow_format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        PlatformConfig(),
        TimestampRange(0, 1));

    auto snda = SOMASparseNDArray::open(
        uri, OpenMode::write, ctx, TimestampRange(0, 2));

    int32_t val = 100;
    snda->set_metadata("md", TILEDB_INT32, 1, &val);
    snda->close();

    // Read metadata
    snda->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(snda->metadata_num() == 3);
    REQUIRE(snda->has_metadata("soma_object_type"));
    REQUIRE(snda->has_metadata("soma_encoding_version"));
    REQUIRE(snda->has_metadata("md"));
    auto mdval = snda->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    snda->close();

    // md should not be available at (0, 1)
    snda->open(OpenMode::read, TimestampRange(0, 1));
    REQUIRE(snda->metadata_num() == 2);
    REQUIRE(snda->has_metadata("soma_object_type"));
    REQUIRE(snda->has_metadata("soma_encoding_version"));
    REQUIRE(!snda->has_metadata("md"));
    snda->close();

    // Metadata should also be retrievable in write mode
    snda->open(OpenMode::write);
    REQUIRE(snda->metadata_num() == 3);
    REQUIRE(snda->has_metadata("soma_object_type"));
    REQUIRE(snda->has_metadata("soma_encoding_version"));
    REQUIRE(snda->has_metadata("md"));
    mdval = snda->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    snda->delete_metadata("md");
    REQUIRE(!snda->has_metadata("md"));
    REQUIRE(snda->metadata_num() == 2);
    mdval = snda->get_metadata("md");
    REQUIRE(!mdval.has_value());
    snda->close();

    // Confirm delete in read mode
    snda->open(OpenMode::read);
    REQUIRE(!snda->has_metadata("md"));
    REQUIRE(snda->metadata_num() == 2);
}

TEST_CASE(
    "SOMASparseNDArray: can_tiledbsoma_upgrade_shape", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-sparse-ndarray-upgrade-shape";

    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        attr_tiledb_datatype);

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(
        uri,
        attr_arrow_format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx);

    auto snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
    REQUIRE(snda->has_current_domain());

    auto dom = snda->soma_domain_slot<int64_t>(dim_name);
    REQUIRE(dom.first == 0);
    REQUIRE(dom.second == dim_max);

    std::vector<int64_t> newshape_wrong_dims({dim_max, 12});
    std::vector<int64_t> newshape_too_big({dim_max + 10});
    std::vector<int64_t> newshape_good({40});

    auto check = snda->can_upgrade_shape(newshape_wrong_dims, "testing");
    REQUIRE(check.first == false);
    REQUIRE(
        check.second ==
        "testing: provided shape has ndim 2, while the "
        "array has 1");

    check = snda->can_upgrade_shape(newshape_too_big, "testing");
    REQUIRE(check.first == false);
    REQUIRE(
        check.second ==
        "testing: array already has a shape: please use resize");

    check = snda->can_upgrade_shape(newshape_good, "testing");
    REQUIRE(check.first == false);
    REQUIRE(
        check.second ==
        "testing: array already has a shape: please use resize");
}

TEST_CASE("SOMASparseNDArray: can_resize", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-sparse-ndarray-resize";

    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        attr_tiledb_datatype);

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(
        uri,
        attr_arrow_format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx);

    auto snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
    REQUIRE(snda->has_current_domain());

    // For new-style arrays, with the current-domain feature:
    // * The shape specified at create becomes the core current domain
    //   o Recall that the core current domain is mutable, up tp <= (max) domain
    // * The core (max) domain is huge
    //   o Recall that the core max domain is immutable
    auto dom = snda->soma_domain_slot<int64_t>(dim_name);
    auto mxd = snda->soma_maxdomain_slot<int64_t>(dim_name);
    REQUIRE(dom != mxd);
    REQUIRE(dom.first == 0);
    REQUIRE(dom.second == dim_max);

    std::vector<int64_t> newshape_wrong_dims({dim_max, 12});
    std::vector<int64_t> newshape_too_small({40});
    std::vector<int64_t> newshape_good({2000});

    auto check = snda->can_resize(newshape_wrong_dims, "testing");
    REQUIRE(check.first == false);
    REQUIRE(
        check.second ==
        "testing: provided shape has ndim 2, while the array has 1");

    check = snda->can_resize(newshape_too_small, "testing");
    REQUIRE(check.first == false);
    REQUIRE(
        check.second ==
        "[testing] index-column name 'soma_dim_0': new upper 39 < old upper "
        "999 (downsize is unsupported)");

    check = snda->can_resize(newshape_good, "testing");
    REQUIRE(check.first == true);
    REQUIRE(check.second == "");
}
