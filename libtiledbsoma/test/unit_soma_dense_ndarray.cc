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

TEST_CASE("SOMADenseNDArray: basic", "[SOMADenseNDArray]") {
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
              .string_hi = "N/A",
              .use_current_domain = use_current_domain}});

        auto index_columns = helper::create_column_index_info(dim_infos);

        if (use_current_domain) {
            // Setting a current domain on a TileDB dense array is not (yet)
            // supported
            // https://github.com/single-cell-data/TileDB-SOMA/issues/2955
            REQUIRE_THROWS(SOMADenseNDArray::create(
                uri,
                dim_arrow_format,
                ArrowTable(
                    std::move(index_columns.first),
                    std::move(index_columns.second)),
                ctx,
                PlatformConfig(),
                TimestampRange(0, 2)));
        } else {
            SOMADenseNDArray::create(
                uri,
                attr_arrow_format,
                ArrowTable(
                    std::move(index_columns.first),
                    std::move(index_columns.second)),
                ctx,
                PlatformConfig(),
                TimestampRange(0, 2));

            REQUIRE(SOMADenseNDArray::exists(uri, ctx));
            REQUIRE(!SOMADataFrame::exists(uri, ctx));
            REQUIRE(!SOMASparseNDArray::exists(uri, ctx));

            auto dnda = SOMADenseNDArray::open(uri, OpenMode::read, ctx);
            REQUIRE(dnda->uri() == uri);
            REQUIRE(dnda->ctx() == ctx);
            REQUIRE(dnda->type() == "SOMADenseNDArray");
            REQUIRE(dnda->is_sparse() == false);
            REQUIRE(dnda->soma_data_type() == attr_arrow_format);
            auto schema = dnda->tiledb_schema();
            REQUIRE(schema->has_attribute("soma_data"));
            REQUIRE(schema->array_type() == TILEDB_DENSE);
            REQUIRE(schema->domain().has_dimension(dim_name));
            REQUIRE(dnda->ndim() == 1);

            // TODO: Once we have support for current domain in dense arrays
            // https://github.com/single-cell-data/TileDB-SOMA/issues/2955
            // if (use_current_domain) {
            //    REQUIRE(dnda->shape() == std::vector<int64_t>{dim_max +
            //    1});
            //} else {
            //    REQUIRE(
            //        dnda->maxshape() == std::vector<int64_t>{dim_max +
            //        1});
            //}

            REQUIRE(dnda->maxshape() == std::vector<int64_t>{shape});

            dnda->close();

            std::vector<int64_t> d0{1, 10};
            std::vector<int32_t> a0(10, 1);

            dnda->open(OpenMode::write);
            dnda->set_column_data("soma_data", a0.size(), a0.data());
            dnda->set_column_data(dim_name, d0.size(), d0.data());
            dnda->write();
            dnda->close();

            dnda->open(OpenMode::read);
            while (auto batch = dnda->read_next()) {
                auto arrbuf = batch.value();
                auto a0span = arrbuf->at("soma_data")->data<int32_t>();
                REQUIRE(
                    a0 == std::vector<int32_t>(a0span.begin(), a0span.end()));
            }
            dnda->close();

            dnda = SOMADenseNDArray::open(uri, OpenMode::read, ctx);
            auto new_shape = std::vector<int64_t>({shape * 2});
            // * When use_current_domain is false: can't resize what has not
            //   been sized.
            // * When use_current_domain is true: TODO: current domain not
            //   supported for dense arrays yet (see above).
            REQUIRE_THROWS(dnda->resize(new_shape, "testing"));
            dnda->close();
        }
    }
}

TEST_CASE("SOMADenseNDArray: platform_config", "[SOMADenseNDArray]") {
    int64_t dim_max = 999;
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();
        std::string uri = "mem://unit-test-dense-ndarray-platform-config";
        std::string dim_name = "soma_dim_0";
        tiledb_datatype_t tiledb_datatype = TILEDB_INT64;
        std::string arrow_format = ArrowAdapter::tdb_to_arrow_type(
            tiledb_datatype);

        PlatformConfig platform_config;
        platform_config.dense_nd_array_dim_zstd_level = 6;

        std::vector<helper::DimInfo> dim_infos(
            {{.name = dim_name,
              .tiledb_datatype = tiledb_datatype,
              .dim_max = dim_max,
              .string_lo = "N/A",
              .string_hi = "N/A",
              .use_current_domain = use_current_domain}});

        auto index_columns = helper::create_column_index_info(dim_infos);

        if (use_current_domain) {
            // Setting a current domain on a TileDB dense array is not (yet)
            // supported
            // https://github.com/single-cell-data/TileDB-SOMA/issues/2955
            REQUIRE_THROWS(SOMADenseNDArray::create(
                uri,
                arrow_format,
                ArrowTable(
                    std::move(index_columns.first),
                    std::move(index_columns.second)),
                ctx,
                platform_config));

        } else {
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
            REQUIRE(
                dim_filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL) == 6);

            dnda->close();
        }
    }
}

TEST_CASE("SOMADenseNDArray: metadata", "[SOMADenseNDArray]") {
    int64_t dim_max = 999;
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();
        std::string uri = "mem://unit-test-dense-ndarray";
        std::string dim_name = "soma_dim_0";
        tiledb_datatype_t tiledb_datatype = TILEDB_INT64;
        std::string arrow_format = ArrowAdapter::tdb_to_arrow_type(
            tiledb_datatype);

        std::vector<helper::DimInfo> dim_infos(
            {{.name = dim_name,
              .tiledb_datatype = tiledb_datatype,
              .dim_max = dim_max,
              .string_lo = "N/A",
              .string_hi = "N/A",
              .use_current_domain = use_current_domain}});

        auto index_columns = helper::create_column_index_info(dim_infos);

        SOMASparseNDArray::create(
            uri,
            arrow_format,
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            PlatformConfig(),
            TimestampRange(0, 2));

        auto dnda = SOMADenseNDArray::open(
            uri,
            OpenMode::write,
            ctx,
            {},
            ResultOrder::automatic,
            std::pair<uint64_t, uint64_t>(1, 1));

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
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
        dnda->close();

        // md should not be available at (2, 2)
        dnda->open(OpenMode::read, TimestampRange(2, 2));
        REQUIRE(dnda->metadata_num() == 2);
        REQUIRE(dnda->has_metadata("soma_object_type"));
        REQUIRE(dnda->has_metadata("soma_encoding_version"));
        REQUIRE(!dnda->has_metadata("md"));
        dnda->close();

        // Metadata should also be retrievable in write mode
        dnda->open(OpenMode::write, TimestampRange(0, 2));
        REQUIRE(dnda->metadata_num() == 3);
        REQUIRE(dnda->has_metadata("soma_object_type"));
        REQUIRE(dnda->has_metadata("soma_encoding_version"));
        REQUIRE(dnda->has_metadata("md"));
        mdval = dnda->get_metadata("md");
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

        // Delete and have it reflected when reading metadata while in write
        // mode
        dnda->delete_metadata("md");
        mdval = dnda->get_metadata("md");
        REQUIRE(!mdval.has_value());
        dnda->close();

        // Confirm delete in read mode
        dnda->open(OpenMode::read, TimestampRange(0, 2));
        REQUIRE(!dnda->has_metadata("md"));
        REQUIRE(dnda->metadata_num() == 2);
    }
}
