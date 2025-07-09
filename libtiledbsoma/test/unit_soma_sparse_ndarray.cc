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

#include <stdexcept>
#include <utility>
#include <variant>

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
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(attr_tiledb_datatype);

    REQUIRE(!SOMASparseNDArray::exists(uri, ctx));

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(uri, attr_arrow_format, index_columns, ctx, PlatformConfig(), TimestampRange(0, 2));

    REQUIRE(SOMASparseNDArray::exists(uri, ctx));
    REQUIRE(!SOMADataFrame::exists(uri, ctx));
    REQUIRE(!SOMADenseNDArray::exists(uri, ctx));

    auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_read, ctx);
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
    {
        snda->open(OpenMode::soma_read);
        auto mq = snda->create_managed_query();
        mq.setup_write_column(dim_name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column(attr_name, a0.size(), a0.data(), (uint64_t*)nullptr);
        REQUIRE_THROWS(mq.submit_write());
        snda->close();
    }

    {
        snda->open(OpenMode::soma_write);
        auto mq = snda->create_managed_query();
        mq.setup_write_column(dim_name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column(attr_name, a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        snda->close();
    }

    {
        snda->open(OpenMode::soma_read);
        auto mq = snda->create_managed_query();
        while (auto batch = mq.read_next()) {
            auto arrbuf = batch.value();
            auto d0span = arrbuf->at(dim_name)->data<int64_t>();
            auto a0span = arrbuf->at(attr_name)->data<int32_t>();
            REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
            REQUIRE(a0 == std::vector<int32_t>(a0span.begin(), a0span.end()));
        }
        snda->close();
    }

    std::vector<int64_t> d0b({dim_max, dim_max + 1});
    std::vector<int64_t> a0b({30, 40});

    // Try out-of-bounds write before resize.
    // * Without current domain support: this should throw since it's
    //   outside the (immutable) doqain.
    // * With current domain support: this should throw since it's outside
    // the (mutable) current domain.
    {
        snda->open(OpenMode::soma_write);
        auto mq = snda->create_managed_query();
        mq.setup_write_column(dim_name, d0b.size(), d0b.data(), (uint64_t*)nullptr);
        mq.setup_write_column(attr_name, a0b.size(), a0b.data(), (uint64_t*)nullptr);
        REQUIRE_THROWS(mq.submit_write());
        snda->close();
    }

    auto new_shape = std::vector<int64_t>({shape * 2});

    snda->open(OpenMode::soma_write);
    // Should throw since this already has a shape (core current
    // domain).
    REQUIRE_THROWS(snda->upgrade_shape(new_shape, "testing"));
    snda->resize(new_shape, "testing");
    snda->close();

    // Try out-of-bounds write after resize.
    {
        snda->open(OpenMode::soma_write);
        auto mq = snda->create_managed_query();
        mq.setup_write_column(dim_name, d0b.size(), d0b.data(), (uint64_t*)nullptr);
        mq.setup_write_column(attr_name, a0b.size(), a0b.data(), (uint64_t*)nullptr);
        // Implicitly checking for no throw
        mq.submit_write();
        snda->close();
    }
    snda->open(OpenMode::soma_read);
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
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(attr_tiledb_datatype);

    PlatformConfig platform_config;
    platform_config.sparse_nd_array_dim_zstd_level = 6;

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(uri, attr_arrow_format, index_columns, ctx, platform_config);

    auto soma_dataframe = SOMASparseNDArray::open(uri, OpenMode::soma_read, ctx);
    auto dim_filter = soma_dataframe->tiledb_schema()->domain().dimension(dim_name).filter_list().filter(0);
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
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(attr_tiledb_datatype);

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(uri, attr_arrow_format, index_columns, ctx, PlatformConfig(), TimestampRange(0, 1));

    auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx, TimestampRange(0, 2));

    int32_t val = 100;
    snda->set_metadata("md", TILEDB_INT32, 1, &val);
    snda->close();

    // Read metadata
    snda->open(OpenMode::soma_read, TimestampRange(0, 2));
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
    snda->open(OpenMode::soma_read, TimestampRange(0, 1));
    REQUIRE(snda->metadata_num() == 2);
    REQUIRE(snda->has_metadata("soma_object_type"));
    REQUIRE(snda->has_metadata("soma_encoding_version"));
    REQUIRE(!snda->has_metadata("md"));
    snda->close();

    // Metadata should also be retrievable in write mode
    snda->open(OpenMode::soma_write);
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
    snda->open(OpenMode::soma_read);
    REQUIRE(!snda->has_metadata("md"));
    REQUIRE(snda->metadata_num() == 2);
}

TEST_CASE("SOMASparseNDArray: can_tiledbsoma_upgrade_shape", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-sparse-ndarray-upgrade-shape";

    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(attr_tiledb_datatype);

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(uri, attr_arrow_format, index_columns, ctx);

    auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx);
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
    REQUIRE(check.second == "testing: array already has a shape: please use resize");

    check = snda->can_upgrade_shape(newshape_good, "testing");
    REQUIRE(check.first == false);
    REQUIRE(check.second == "testing: array already has a shape: please use resize");
}

TEST_CASE("SOMASparseNDArray: can_resize", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-sparse-ndarray-resize";

    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(attr_tiledb_datatype);

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(uri, attr_arrow_format, index_columns, ctx);

    auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx);
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
    REQUIRE(check.second == "testing: provided shape has ndim 2, while the array has 1");

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

TEST_CASE("SOMASparseNDArray: nnz", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-sparse-ndarray-resize";

    std::string dim_name = "soma_dim_0";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(attr_tiledb_datatype);

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(uri, attr_arrow_format, index_columns, ctx);

    // Create vectors of data for writing.
    std::vector<int64_t> d0 = {1, 2, 3};
    std::vector<int32_t> a0 = {0, 0, 0};

    {
        auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx);
        auto mq = snda->create_managed_query();
        mq.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column("soma_data", a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        snda->close();

        REQUIRE(snda->nnz() == 3);
    }
    // Write to non overlapping fragment
    d0 = {4, 5, 6};
    a0 = {1, 1, 1};

    {
        auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx);
        auto mq = snda->create_managed_query();
        mq.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column("soma_data", a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        snda->close();

        REQUIRE(snda->nnz() == 6);
    }

    // Write to partially overlapping fragment
    d0 = {6, 7, 8};
    a0 = {2, 2, 2};

    {
        auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx);
        auto mq = snda->create_managed_query();
        mq.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column("soma_data", a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        snda->close();

        REQUIRE(snda->nnz() == 8);
    }

    // Write to overlapping fragment
    d0 = {1, 4, 8};
    a0 = {3, 3, 3};
    {
        auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx);
        auto mq = snda->create_managed_query();
        mq.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column("soma_data", a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        snda->close();

        REQUIRE(snda->nnz() == 8);
    }

    // Write separate overlapping fragments
    d0 = {10, 11, 12};
    a0 = {4, 4, 4};

    {
        auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx);
        auto mq1 = snda->create_managed_query();
        mq1.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq1.setup_write_column("soma_data", a0.size(), a0.data(), (uint64_t*)nullptr);
        mq1.submit_write();

        mq1.reset();

        d0 = {12, 13, 14};
        a0 = {5, 5, 5};
        auto mq2 = snda->create_managed_query();
        mq2.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq2.setup_write_column("soma_data", a0.size(), a0.data(), (uint64_t*)nullptr);
        mq2.submit_write();
        snda->close();

        REQUIRE(snda->nnz() == 13);
    }
}

TEST_CASE("SOMASparseNDArray: resize with timestamp", "[SOMASparseNDArray]") {
    // Core uses domain & current domain like (0, 999); SOMA uses shape like
    // 1000. We want to carefully and explicitly test here that there aren't any
    // off-by-one errors.
    int64_t dim_max = 999;
    int64_t shape = 1000;

    auto orig_shape = std::vector<int64_t>({shape});
    auto new_shape = std::vector<int64_t>({shape * 2});

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-sparse-ndarray-basic";
    std::string dim_name = "soma_dim_0";
    std::string attr_name = "soma_data";
    tiledb_datatype_t dim_tiledb_datatype = TILEDB_INT64;
    tiledb_datatype_t attr_tiledb_datatype = TILEDB_INT32;
    std::string dim_arrow_format = ArrowAdapter::tdb_to_arrow_type(dim_tiledb_datatype);
    std::string attr_arrow_format = ArrowAdapter::tdb_to_arrow_type(attr_tiledb_datatype);

    REQUIRE(!SOMASparseNDArray::exists(uri, ctx));

    std::vector<helper::DimInfo> dim_infos(
        {{.name = dim_name,
          .tiledb_datatype = dim_tiledb_datatype,
          .dim_max = dim_max,
          .string_lo = "N/A",
          .string_hi = "N/A"}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(uri, attr_arrow_format, index_columns, ctx, PlatformConfig(), TimestampRange(0, 1));

    index_columns.first->release(index_columns.first.get());
    index_columns.second->release(index_columns.second.get());

    auto snda = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx, TimestampRange(0, 2));
    snda->resize(new_shape, "testing");
    snda->close();

    snda->open(OpenMode::soma_read, TimestampRange(0, 1));
    REQUIRE(snda->shape() == orig_shape);
    snda->close();

    snda->open(OpenMode::soma_read, TimestampRange(0, 2));
    REQUIRE(snda->shape() == new_shape);
    snda->close();
}

TEST_CASE("SOMASparseNDArray: delete cells", "[SOMASparseNDArray][delete]") {
    // Create a TileDB sparse array with 2 integer dimensions.
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://test-query-condition-sparse";
    auto index_columns = helper::create_column_index_info(
        {helper::DimInfo(
             {.name = "soma_dim_0",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = 3,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "soma_dim_1",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = 2,
              .string_lo = "N/A",
              .string_hi = "N/A"})});

    SOMASparseNDArray::create(uri, "i", index_columns, ctx);

    {
        INFO("Write data to array.");
        std::vector<int64_t> coords_dim_0{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
        std::vector<int64_t> coords_dim_1{0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
        std::vector<int32_t> data(12);
        std::iota(data.begin(), data.end(), 1);

        Array array{*ctx->tiledb_ctx(), uri, TILEDB_WRITE};
        Query query{*ctx->tiledb_ctx(), array};
        query.set_layout(TILEDB_GLOBAL_ORDER);
        query.set_data_buffer("soma_dim_0", coords_dim_0);
        query.set_data_buffer("soma_dim_1", coords_dim_1);
        query.set_data_buffer("soma_data", data);
        query.submit();
        query.finalize();
        array.close();
    }

    // Create variable for tests. Using sections will rerun the this test from beginning to end for each section.
    std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{};
    int64_t expected_result_num{0};
    std::vector<int32_t> expected_data{};
    std::vector<int64_t> expected_dim_0{};
    std::vector<int64_t> expected_dim_1{};

    auto check_delete = [&](const std::string& log_note) {
        INFO(log_note);
        expected_data.resize(12, 0);
        expected_dim_0.resize(12, 0);
        expected_dim_1.resize(12, 0);

        {
            INFO("Delete cells from the sparse array.");
            auto sparse_array = SOMASparseNDArray::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
            sparse_array->delete_cells(delete_coords);
            sparse_array->close();
        }

        std::vector<int32_t> actual_data(12);
        std::vector<int64_t> actual_dim_0(12);
        std::vector<int64_t> actual_dim_1(12);
        Array array{*ctx->tiledb_ctx(), uri, TILEDB_READ};
        Query query{*ctx->tiledb_ctx(), array};
        Subarray subarray(*ctx->tiledb_ctx(), array);
        subarray.add_range<int64_t>(0, 0, 3).add_range<int64_t>(1, 0, 2);
        query.set_layout(TILEDB_GLOBAL_ORDER);
        query.set_subarray(subarray);
        query.set_data_buffer("soma_dim_0", actual_dim_0);
        query.set_data_buffer("soma_dim_1", actual_dim_1);
        query.set_data_buffer("soma_data", actual_data);
        query.submit();
        query.finalize();
        array.close();

        auto actual_result_num = static_cast<int64_t>(query.result_buffer_elements()["soma_data"].second);
        CHECK(actual_result_num == expected_result_num);
        CHECK_THAT(actual_dim_0, Catch::Matchers::Equals(expected_dim_0));
        CHECK_THAT(actual_dim_1, Catch::Matchers::Equals(expected_dim_1));
        CHECK_THAT(actual_data, Catch::Matchers::Equals(actual_data));
    };

    SECTION("Delete all using ranges") {
        expected_result_num = 0;
        delete_coords.assign({std::pair<int64_t, int64_t>(0, 3), std::pair<int64_t, int64_t>(0, 2)});
        check_delete("Delete all using ranges");
    }
    SECTION("Delete all using row ranges") {
        expected_result_num = 0;
        delete_coords.assign({std::pair<int64_t, int64_t>(0, 3), std::pair<int64_t, int64_t>(0, 2)});
        check_delete("Delete all using row range");
    }
    SECTION("Delete all using column range") {
        expected_result_num = 0;
        delete_coords.assign({std::monostate(), std::pair<int64_t, int64_t>(0, 2)});
        check_delete("Delete all using column range");
    }
    SECTION("Delete all using coordinates") {
        expected_result_num = 0;
        delete_coords.assign({std::vector<int64_t>({0, 1, 2, 3}), std::vector<int64_t>({0, 1, 2})});
        check_delete("Delete all using coordinates");
    }
    SECTION("Delete 1 row using range") {
        expected_result_num = 9;
        delete_coords.assign({std::vector<int64_t>({1, 1})});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using range");
    }
    SECTION("Delete 1 row using coordinate") {
        expected_result_num = 9;
        delete_coords.assign({std::vector<int64_t>({1})});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using coordinate");
    }
    SECTION("Delete 1 row using row range, empty coord (select all)") {
        expected_result_num = 9;
        delete_coords.assign({std::pair<int64_t, int64_t>({1, 1}), std::monostate()});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using row range, empty coord");
    }
    SECTION("Delete 1 row using duplicate coordinates") {
        expected_result_num = 9;
        delete_coords.assign({std::vector<int64_t>({1, 1, 1})});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using duplicate coordinates");
    }
    SECTION("Delete 1 row using a coordinate") {
        expected_result_num = 9;
        delete_coords.assign({std::vector<int64_t>({1})});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using a coordinate");
    }
    SECTION("Delete multiple rows with coordinates (unordered)") {
        expected_result_num = 3;
        delete_coords.assign({{std::vector<int64_t>({3, 0, 1})}});
        expected_data.assign({7, 8, 9});
        expected_dim_0.assign({2, 2, 2});
        expected_dim_1.assign({0, 1, 2});
        check_delete("Delete multiple rows with coordinates (unordered)");
    }
    SECTION("Delete range on row, range on column") {
        expected_result_num = 8;
        delete_coords.assign({std::pair<int64_t, int64_t>({0, 1}), std::pair<int64_t, int64_t>({1, 2})});
        expected_data.assign({1, 4, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 1, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 0, 0, 1, 2, 0, 1, 2});
        check_delete("Delete range on row, coord on column");
    }
    SECTION("Delete range on row, coords on column") {
        expected_result_num = 9;
        delete_coords.assign({std::pair<int64_t, int64_t>({0, 2}), std::vector<int64_t>({1})});
        expected_data.assign({1, 3, 4, 6, 7, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 1, 1, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 2, 0, 2, 0, 2, 0, 1, 2});
        check_delete("Delete range on row, coords on column");
    }
    SECTION("Delete coords on row, range on column") {
        expected_result_num = 8;
        delete_coords.assign({std::vector<int64_t>({1, 3}), std::pair<int64_t, int64_t>({0, 1})});
        expected_data.assign({1, 2, 3, 6, 7, 8, 9, 12});
        expected_dim_0.assign({0, 0, 0, 1, 2, 2, 2, 3});
        expected_dim_1.assign({0, 1, 2, 2, 0, 1, 2, 2});
        check_delete("Delete coords on row, range on column");
    }
    SECTION("Delete coords on row, coords on column") {
        expected_result_num = 6;
        delete_coords.assign({std::vector<int64_t>({3, 0, 2}), std::vector<int64_t>({0, 2})});
        expected_data.assign({2, 4, 5, 6, 8, 11});
        expected_dim_0.assign({0, 1, 1, 1, 2, 3});
        expected_dim_1.assign({1, 0, 1, 2, 1, 1});
        check_delete("Delete coords on row, coords on column");
    }
}

TEST_CASE("SOMASparseNDArray: check delete cell exceptions", "[SOMASparseNDArray][delete]") {
    // Create a TileDB sparse array with 3 integer dimensions.
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://test-query-condition-sparse";
    auto index_columns = helper::create_column_index_info(
        {helper::DimInfo(
             {.name = "soma_dim_0",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = 3,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "soma_dim_1",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = 4,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "soma_dim_2",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = 11,
              .string_lo = "N/A",
              .string_hi = "N/A"})});

    SOMASparseNDArray::create(uri, "i", index_columns, ctx);

    auto sparse_array = SOMASparseNDArray::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
    {
        INFO("Check throws: no coordinates.");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::invalid_argument);
    }
    {
        INFO("Check throws: full range less than current domain (dim=0).");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::pair<int64_t, int64_t>({-10, -1})};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::out_of_range);
    }
    {
        INFO("Check throws: full range less than current domain (dim=1).");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::monostate(), std::pair<int64_t, int64_t>(-10, -1)};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::out_of_range);
    }
    {
        INFO("Check throws: range starts at max value + 1 (dim=0).");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::pair<int64_t, int64_t>(4, 10)};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::out_of_range);
    }
    {
        INFO("Check throws: range starts at max value + 1 (dim=1).");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::pair<int64_t, int64_t>(5, 10)};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::out_of_range);
    }
    {
        INFO("Check throws: invalid range (no values)");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::monostate()};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::invalid_argument);
    }
    {
        INFO("Check throws: invalid range (unordered)");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::pair<int64_t, int64_t>(3, 1)};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::invalid_argument);
    }
    {
        INFO("Check throws: invalid range (unordered) that is also out of bounds");
        // Note: the invalid argument error should take precedent over the fact the range is out-of-bounds.
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::pair<int64_t, int64_t>(20, 11), std::monostate()};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::invalid_argument);
    }
    {
        INFO("Check throws: coordinate out-of-bounds (dim=0) ");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::vector<int64_t>({1, 100, 2}), std::vector<int64_t>({1, 3})};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::out_of_range);
    }
    {
        INFO("Check throws: coordinate out-of-bounds (dim=2)");
        std::vector<std::variant<std::monostate, std::pair<int64_t, int64_t>, std::vector<int64_t>>> delete_coords{
            std::pair<int64_t, int64_t>(1, 1), std::monostate(), std::vector<int64_t>({-1, 1, 4, 2, 11})};
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_coords), std::out_of_range);
    }
    sparse_array->close();
}
