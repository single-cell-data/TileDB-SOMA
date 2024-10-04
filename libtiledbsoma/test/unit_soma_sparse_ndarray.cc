/**
 * @file   unit_soma_sparse_ndarray.cc
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
 * This file manages unit tests for the SOMASparseNDArray class
 */

#include "common.h"

TEST_CASE("SOMASparseNDArray: basic", "[SOMASparseNDArray]") {
    // Core uses domain & current domain like (0, 999); SOMA uses shape like
    // 1000. We want to carefully and explicitly test here that there aren't any
    // off-by-one errors.
    int64_t dim_max = 999;
    int64_t shape = 1000;

    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
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
              .string_hi = "N/A",
              .use_current_domain = use_current_domain}});

        auto index_columns = helper::create_column_index_info(dim_infos);

        SOMASparseNDArray::create(
            uri,
            attr_arrow_format,
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
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
        if (!use_current_domain) {
            REQUIRE(snda->maxshape() == expect);
        }

        snda->close();

        std::vector<int64_t> d0(10);
        for (int j = 0; j < 10; j++)
            d0[j] = j;
        std::vector<int32_t> a0(10, 1);

        snda->open(OpenMode::write);
        snda->set_column_data(dim_name, d0.size(), d0.data());
        snda->set_column_data(attr_name, a0.size(), a0.data());
        snda->write();
        snda->close();

        snda->open(OpenMode::read);
        while (auto batch = snda->read_next()) {
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
        snda->set_column_data(dim_name, d0b.size(), d0b.data());
        snda->set_column_data(attr_name, a0b.size(), a0b.data());
        REQUIRE_THROWS(snda->write());
        snda->close();

        if (!use_current_domain) {
            auto new_shape = std::vector<int64_t>({shape});

            snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
            // Without current-domain support: this should throw since
            // one cannot resize what has not been sized.
            REQUIRE(!snda->has_current_domain());
            REQUIRE_THROWS(snda->resize(new_shape, "testing"));
            // Now set the shape
            snda->upgrade_shape(new_shape, "testing");
            snda->close();

            snda->open(OpenMode::read);
            REQUIRE(snda->has_current_domain());
            snda->close();

            snda->open(OpenMode::write);
            REQUIRE(snda->has_current_domain());
            // Should not fail since we're setting it to what it already is.
            snda->resize(new_shape, "testing");
            snda->close();

            snda = SOMASparseNDArray::open(uri, OpenMode::read, ctx);
            REQUIRE(snda->shape() == new_shape);
            snda->close();

        } else {
            auto new_shape = std::vector<int64_t>({shape * 2});

            snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
            // Should throw since this already has a shape (core current
            // domain).
            REQUIRE_THROWS(snda->upgrade_shape(new_shape, "testing"));
            snda->resize(new_shape, "testing");
            snda->close();

            // Try out-of-bounds write after resize.
            snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
            snda->set_column_data(dim_name, d0b.size(), d0b.data());
            snda->set_column_data(attr_name, a0b.size(), a0b.data());
            // Implicitly checking for no throw
            snda->write();
            snda->close();

            snda->open(OpenMode::read);
            REQUIRE(snda->shape() == new_shape);
            snda->close();
        }
    }
}

TEST_CASE("SOMASparseNDArray: platform_config", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
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
              .string_hi = "N/A",
              .use_current_domain = use_current_domain}});

        auto index_columns = helper::create_column_index_info(dim_infos);

        SOMASparseNDArray::create(
            uri,
            attr_arrow_format,
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
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
}

TEST_CASE("SOMASparseNDArray: metadata", "[SOMASparseNDArray]") {
    int64_t dim_max = 999;
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
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
              .string_hi = "N/A",
              .use_current_domain = use_current_domain}});

        auto index_columns = helper::create_column_index_info(dim_infos);

        SOMASparseNDArray::create(
            uri,
            attr_arrow_format,
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            PlatformConfig(),
            TimestampRange(0, 2));

        auto snda = SOMASparseNDArray::open(
            uri,
            OpenMode::write,
            ctx,
            {},
            ResultOrder::automatic,
            std::pair<uint64_t, uint64_t>(1, 1));

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
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
        snda->close();

        // md should not be available at (2, 2)
        snda->open(OpenMode::read, TimestampRange(2, 2));
        REQUIRE(snda->metadata_num() == 2);
        REQUIRE(snda->has_metadata("soma_object_type"));
        REQUIRE(snda->has_metadata("soma_encoding_version"));
        REQUIRE(!snda->has_metadata("md"));
        snda->close();

        // Metadata should also be retrievable in write mode
        snda->open(OpenMode::write, TimestampRange(0, 2));
        REQUIRE(snda->metadata_num() == 3);
        REQUIRE(snda->has_metadata("soma_object_type"));
        REQUIRE(snda->has_metadata("soma_encoding_version"));
        REQUIRE(snda->has_metadata("md"));
        mdval = snda->get_metadata("md");
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

        // Delete and have it reflected when reading metadata while in write
        // mode
        snda->delete_metadata("md");
        mdval = snda->get_metadata("md");
        REQUIRE(!mdval.has_value());
        snda->close();

        // Confirm delete in read mode
        snda->open(OpenMode::read, TimestampRange(0, 2));
        REQUIRE(!snda->has_metadata("md"));
        REQUIRE(snda->metadata_num() == 2);
    }
}
void breakme() {
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
          .string_hi = "N/A",
          .use_current_domain = false}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(
        uri,
        attr_arrow_format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx);

    auto snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
    REQUIRE(snda->has_current_domain() == false);

    // For old-style arrays, from before the current-domain feature:
    // * The shape specified at create becomes the core (max) domain
    //   o Recall that the core domain is immutable
    // * There is no current domain set
    //   o A current domain can be applied to it, up to <= (max) domain
    auto dom = snda->soma_domain_slot<int64_t>(dim_name);
    auto mxd = snda->soma_maxdomain_slot<int64_t>(dim_name);
    REQUIRE(dom == mxd);
    REQUIRE(dom.first == 0);
    REQUIRE(dom.second == dim_max);

    breakme();
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
        "testing for soma_dim_0: new 1009 < maxshape "
        "1000");

    check = snda->can_upgrade_shape(newshape_good, "testing");
    REQUIRE(check.first == true);
    REQUIRE(check.second == "");
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
          .string_hi = "N/A",
          .use_current_domain = true}});

    auto index_columns = helper::create_column_index_info(dim_infos);

    SOMASparseNDArray::create(
        uri,
        attr_arrow_format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx);

    auto snda = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
    REQUIRE(snda->has_current_domain() == true);

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
        check.second == "testing for soma_dim_0: new 40 < existing shape 1000");

    check = snda->can_resize(newshape_good, "testing");
    REQUIRE(check.first == true);
    REQUIRE(check.second == "");
}
