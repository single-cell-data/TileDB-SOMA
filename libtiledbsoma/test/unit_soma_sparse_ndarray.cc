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
    CHECK(snda->uri() == uri);
    CHECK(snda->ctx() == ctx);
    CHECK(snda->ndim() == 1);
    CHECK(snda->type() == "SOMASparseNDArray");
    CHECK(snda->is_sparse() == true);
    CHECK(snda->soma_data_type() == attr_arrow_format);

    auto schema = snda->tiledb_schema();
    CHECK(schema->array_type() == TILEDB_SPARSE);
    // TODO: Check capacity, tile/cell order, etc.

    CHECK(schema->attribute_num() == 1);
    REQUIRE(schema->has_attribute(attr_name));
    auto attr = schema->attribute(attr_name);
    CHECK(attr.type() == TILEDB_INT32);

    auto domain = schema->domain();
    CHECK(domain.ndim() == 1);
    REQUIRE(domain.has_dimension(dim_name));
    auto dim = domain.dimension(dim_name);
    REQUIRE(dim.type() == TILEDB_INT64);
    auto dim_domain = dim.domain<int64_t>();
    CHECK(dim_domain.first == 0);
    CHECK(dim_domain.second == 2147483646);
    CHECK(dim.tile_extent<int64_t>() == 1);

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
        REQUIRE_THROWS(mq.setup_write_column(dim_name, d0.size(), d0.data(), (uint64_t*)nullptr));
        REQUIRE_THROWS(mq.setup_write_column(attr_name, a0.size(), a0.data(), (uint64_t*)nullptr));
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

    auto shape = snda->shape();
    REQUIRE(shape[0] - 1 == dim_max);

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
    auto shape = snda->shape();
    auto max_shape = snda->maxshape();
    REQUIRE(shape != max_shape);
    REQUIRE(shape[0] - 1 == dim_max);

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
        Array array{*ctx->tiledb_ctx(), uri, TILEDB_WRITE};
        Query query{*ctx->tiledb_ctx(), array};
        query.set_layout(TILEDB_GLOBAL_ORDER);

        {
            std::vector<int64_t> coords_dim_0{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
            std::vector<int64_t> coords_dim_1{0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
            std::vector<int32_t> data(12);
            std::iota(data.begin(), data.end(), 1);

            query.set_data_buffer("soma_dim_0", coords_dim_0);
            query.set_data_buffer("soma_dim_1", coords_dim_1);
            query.set_data_buffer("soma_data", data);
            query.submit();
        }

        query.finalize();
        array.close();
    }

    // Create variable for tests. Using sections will rerun the this test from beginning to end for each section.
    std::vector<SOMAColumnSelection<int64_t>> delete_coords{};
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
            auto delete_filter = sparse_array->create_coordinate_value_filter();
            for (size_t index = 0; index < delete_coords.size(); ++index) {
                delete_filter.add_column_selection<int64_t>(index, delete_coords[index]);
            }
            sparse_array->delete_cells(delete_filter);
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
        CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
    };

    SECTION("Delete all using ranges") {
        expected_result_num = 0;
        delete_coords.assign({SliceSelection<int64_t>(0, 3), SliceSelection<int64_t>(0, 2)});
        check_delete("Delete all using ranges");
    }
    SECTION("Delete all using row ranges") {
        expected_result_num = 0;
        delete_coords.assign({SliceSelection<int64_t>(0, 3), SliceSelection<int64_t>(0, 2)});
        check_delete("Delete all using row range");
    }
    SECTION("Delete all using column range") {
        expected_result_num = 0;
        delete_coords.assign({std::monostate(), SliceSelection<int64_t>(0, 2)});
        check_delete("Delete all using column range");
    }
    SECTION("Delete all using coordinates") {
        expected_result_num = 0;
        std::vector<int64_t> points1{0, 1, 2, 3};
        std::vector<int64_t> points2{0, 1, 2};
        delete_coords.assign({PointSelection<int64_t>(points1), PointSelection<int64_t>(points2)});
        check_delete("Delete all using coordinates");
    }
    SECTION("Delete 1 row using range") {
        expected_result_num = 9;
        std::vector<int64_t> points1{1, 1};
        delete_coords.assign({PointSelection<int64_t>(points1)});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using range");
    }
    SECTION("Delete 1 row using coordinate") {
        expected_result_num = 9;
        std::vector<int64_t> points{1};
        delete_coords.assign({PointSelection<int64_t>(points)});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using coordinate");
    }
    SECTION("Delete 1 row using row range, empty coord (select all)") {
        expected_result_num = 9;
        delete_coords.assign({SliceSelection<int64_t>(1, 1), std::monostate()});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using row range, empty coord");
    }
    SECTION("Delete 1 row using duplicate coordinates") {
        expected_result_num = 9;
        std::vector<int64_t> points{1, 1, 1};
        delete_coords.assign({PointSelection<int64_t>(points)});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using duplicate coordinates");
    }
    SECTION("Delete 1 row using a coordinate") {
        expected_result_num = 9;
        std::vector<int64_t> points{1};
        delete_coords.assign({PointSelection<int64_t>(points)});
        expected_data.assign({1, 2, 3, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 0, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 1, 2, 0, 1, 2, 0, 1, 2});
        check_delete("Delete 1 row using a coordinate");
    }
    SECTION("Delete multiple rows with coordinates (unordered)") {
        expected_result_num = 3;
        std::vector<int64_t> points{3, 0, 1};
        delete_coords.assign({PointSelection<int64_t>(points)});
        expected_data.assign({7, 8, 9});
        expected_dim_0.assign({2, 2, 2});
        expected_dim_1.assign({0, 1, 2});
        check_delete("Delete multiple rows with coordinates (unordered)");
    }
    SECTION("Delete range on row, range on column") {
        expected_result_num = 8;
        delete_coords.assign({SliceSelection<int64_t>(0, 1), SliceSelection<int64_t>(1, 2)});
        expected_data.assign({1, 4, 7, 8, 9, 10, 11, 12});
        expected_dim_0.assign({0, 1, 2, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 0, 0, 1, 2, 0, 1, 2});
        check_delete("Delete range on row, coord on column");
    }
    SECTION("Delete range on row, coords on column") {
        expected_result_num = 9;
        std::vector<int64_t> points{1};
        delete_coords.assign({SliceSelection<int64_t>(0, 2), PointSelection<int64_t>(points)});
        expected_data.assign({1, 3, 4, 6, 7, 9, 10, 11, 12});
        expected_dim_0.assign({0, 0, 1, 1, 2, 2, 3, 3, 3});
        expected_dim_1.assign({0, 2, 0, 2, 0, 2, 0, 1, 2});
        check_delete("Delete range on row, coords on column");
    }
    SECTION("Delete coords on row, range on column") {
        expected_result_num = 8;
        std::vector<int64_t> points{1, 3};
        delete_coords.assign({PointSelection<int64_t>(points), SliceSelection<int64_t>(0, 1)});
        expected_data.assign({1, 2, 3, 6, 7, 8, 9, 12});
        expected_dim_0.assign({0, 0, 0, 1, 2, 2, 2, 3});
        expected_dim_1.assign({0, 1, 2, 2, 0, 1, 2, 2});
        check_delete("Delete coords on row, range on column");
    }
    SECTION("Delete coords on row, coords on column") {
        expected_result_num = 6;
        std::vector<int64_t> points1{3, 0, 2};
        std::vector<int64_t> points2{0, 2};
        delete_coords.assign({PointSelection<int64_t>(points1), PointSelection<int64_t>(points2)});
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

    {
        INFO("Check cannot delete in read mode.");
        auto sparse_array = SOMASparseNDArray::open(uri, OpenMode::soma_read, ctx, std::nullopt);
        auto delete_filters = sparse_array->create_coordinate_value_filter();
        delete_filters.add_slice<int64_t>(0, SliceSelection<int64_t>(0, 1));
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_filters), TileDBSOMAError);
        sparse_array->close();
    }

    {
        INFO("Check cannot delete in write mode.");
        auto sparse_array = SOMASparseNDArray::open(uri, OpenMode::soma_write, ctx, std::nullopt);
        auto delete_filters = sparse_array->create_coordinate_value_filter();
        delete_filters.add_slice<int64_t>(0, SliceSelection<int64_t>(0, 1));
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_filters), TileDBSOMAError);
        sparse_array->close();
    }

    auto sparse_array = SOMASparseNDArray::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
    {
        INFO("Check throws: no coordinates.");
        auto delete_filter = sparse_array->create_coordinate_value_filter();
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_filter), std::invalid_argument);
    }
    {
        INFO("Check throws: invalid range (no values)");
        auto delete_filter = sparse_array->create_coordinate_value_filter();
        delete_filter.add_column_selection<int64_t>(0, std::monostate());
        CHECK_THROWS_AS(sparse_array->delete_cells(delete_filter), std::invalid_argument);
    }
    sparse_array->close();
}
