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
    int64_t dim_max = 1000;
    // auto use_current_domain = GENERATE(false, true);
    auto use_current_domain = GENERATE(true);
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
        tiledb_datatype_t tiledb_datatype = TILEDB_INT64;
        std::string arrow_format = ArrowAdapter::tdb_to_arrow_type(
            tiledb_datatype);

        REQUIRE(!SOMASparseNDArray::exists(uri, ctx));

        std::vector<helper::DimInfo> dim_infos(
            {{.name = dim_name,
              .tiledb_datatype = tiledb_datatype,
              .dim_max = dim_max,
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

        REQUIRE(SOMASparseNDArray::exists(uri, ctx));
        REQUIRE(!SOMADataFrame::exists(uri, ctx));
        REQUIRE(!SOMADenseNDArray::exists(uri, ctx));

        auto soma_sparse = SOMASparseNDArray::open(uri, OpenMode::read, ctx);
        REQUIRE(soma_sparse->uri() == uri);
        REQUIRE(soma_sparse->ctx() == ctx);
        REQUIRE(soma_sparse->type() == "SOMASparseNDArray");
        REQUIRE(soma_sparse->is_sparse() == true);
        REQUIRE(soma_sparse->soma_data_type() == arrow_format);
        auto schema = soma_sparse->tiledb_schema();
        REQUIRE(schema->has_attribute(attr_name));
        REQUIRE(schema->array_type() == TILEDB_SPARSE);
        REQUIRE(schema->domain().has_dimension(dim_name));
        REQUIRE(soma_sparse->ndim() == 1);
        REQUIRE(soma_sparse->nnz() == 0);

        auto expect = std::vector<int64_t>({dim_max + 1});
        REQUIRE(soma_sparse->shape() == expect);
        if (!use_current_domain) {
            REQUIRE(soma_sparse->maxshape() == expect);
        }

        soma_sparse->close();

        std::vector<int64_t> d0(10);
        for (int j = 0; j < 10; j++)
            d0[j] = j;
        std::vector<int> a0(10, 1);

        soma_sparse->open(OpenMode::write);
        soma_sparse->set_column_data(dim_name, d0.size(), d0.data());
        soma_sparse->set_column_data(attr_name, a0.size(), a0.data());
        soma_sparse->write();
        soma_sparse->close();

        soma_sparse->open(OpenMode::read);
        while (auto batch = soma_sparse->read_next()) {
            auto arrbuf = batch.value();
            auto d0span = arrbuf->at(dim_name)->data<int64_t>();
            auto a0span = arrbuf->at(attr_name)->data<int>();
            REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
            REQUIRE(a0 == std::vector<int>(a0span.begin(), a0span.end()));
        }

        std::vector<int64_t> d0b({dim_max, dim_max + 1});
        std::vector<int64_t> a0b({30, 40});
        soma_sparse->close();

        // Try out-of-bounds write before resize.
        // * Without current domain support: this should throw since it's
        //   outside the (immutable) doqain.
        // * With current domain support: this should throw since it's outside
        // the (mutable) current domain.
        soma_sparse = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
        soma_sparse->set_column_data(dim_name, d0b.size(), d0b.data());
        soma_sparse->set_column_data(attr_name, a0b.size(), a0b.data());
        REQUIRE_THROWS(soma_sparse->write());
        soma_sparse->close();

        if (!use_current_domain) {
            auto new_shape = std::vector<int64_t>({dim_max});

            soma_sparse = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
            // Without current-domain support: this should throw since
            // one cannot resize what has not been sized.
            REQUIRE(!soma_sparse->has_current_domain());
            REQUIRE_THROWS(soma_sparse->resize(new_shape));
            // Now set the shape
            soma_sparse->upgrade_shape(new_shape);
            REQUIRE(soma_sparse->has_current_domain());
            // Should not fail since we're setting it to what it already is.
            soma_sparse->resize(new_shape);
            soma_sparse->close();

            soma_sparse = SOMASparseNDArray::open(uri, OpenMode::read, ctx);
            REQUIRE(soma_sparse->shape() == new_shape);
            soma_sparse->close();

        } else {
            auto new_shape = std::vector<int64_t>({dim_max * 2});

            soma_sparse = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
            // Should throw since this already has a shape (core current
            // domain).
            REQUIRE_THROWS(soma_sparse->upgrade_shape(new_shape));
            soma_sparse->resize(new_shape);
            soma_sparse->close();

            // Try out-of-bounds write after resize.
            soma_sparse = SOMASparseNDArray::open(uri, OpenMode::write, ctx);
            soma_sparse->set_column_data(dim_name, d0b.size(), d0b.data());
            soma_sparse->set_column_data(attr_name, a0b.size(), a0b.data());
            // Implicitly checking for no throw
            soma_sparse->write();
            soma_sparse->close();
        }
    }
}

TEST_CASE("SOMASparseNDArray: platform_config", "[SOMASparseNDArray]") {
    int64_t dim_max = 1000;
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
        tiledb_datatype_t tiledb_datatype = TILEDB_INT64;
        std::string arrow_format = ArrowAdapter::tdb_to_arrow_type(
            tiledb_datatype);

        PlatformConfig platform_config;
        platform_config.sparse_nd_array_dim_zstd_level = 6;

        std::vector<helper::DimInfo> dim_infos(
            {{.name = dim_name,
              .tiledb_datatype = tiledb_datatype,
              .dim_max = dim_max,
              .use_current_domain = use_current_domain}});

        auto index_columns = helper::create_column_index_info(dim_infos);

        SOMASparseNDArray::create(
            uri,
            arrow_format,
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
    int64_t dim_max = 1000;
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
        tiledb_datatype_t tiledb_datatype = TILEDB_INT64;
        std::string arrow_format = ArrowAdapter::tdb_to_arrow_type(
            tiledb_datatype);

        std::vector<helper::DimInfo> dim_infos(
            {{.name = dim_name,
              .tiledb_datatype = tiledb_datatype,
              .dim_max = dim_max,
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

        auto soma_sparse = SOMASparseNDArray::open(
            uri,
            OpenMode::write,
            ctx,
            {},
            ResultOrder::automatic,
            std::pair<uint64_t, uint64_t>(1, 1));

        int32_t val = 100;
        soma_sparse->set_metadata("md", TILEDB_INT32, 1, &val);
        soma_sparse->close();

        // Read metadata
        soma_sparse->open(OpenMode::read, TimestampRange(0, 2));
        REQUIRE(soma_sparse->metadata_num() == 3);
        REQUIRE(soma_sparse->has_metadata("soma_object_type"));
        REQUIRE(soma_sparse->has_metadata("soma_encoding_version"));
        REQUIRE(soma_sparse->has_metadata("md"));
        auto mdval = soma_sparse->get_metadata("md");
        REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
        REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
        soma_sparse->close();

        // md should not be available at (2, 2)
        soma_sparse->open(OpenMode::read, TimestampRange(2, 2));
        REQUIRE(soma_sparse->metadata_num() == 2);
        REQUIRE(soma_sparse->has_metadata("soma_object_type"));
        REQUIRE(soma_sparse->has_metadata("soma_encoding_version"));
        REQUIRE(!soma_sparse->has_metadata("md"));
        soma_sparse->close();

        // Metadata should also be retrievable in write mode
        soma_sparse->open(OpenMode::write, TimestampRange(0, 2));
        REQUIRE(soma_sparse->metadata_num() == 3);
        REQUIRE(soma_sparse->has_metadata("soma_object_type"));
        REQUIRE(soma_sparse->has_metadata("soma_encoding_version"));
        REQUIRE(soma_sparse->has_metadata("md"));
        mdval = soma_sparse->get_metadata("md");
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

        // Delete and have it reflected when reading metadata while in write
        // mode
        soma_sparse->delete_metadata("md");
        mdval = soma_sparse->get_metadata("md");
        REQUIRE(!mdval.has_value());
        soma_sparse->close();

        // Confirm delete in read mode
        soma_sparse->open(OpenMode::read, TimestampRange(0, 2));
        REQUIRE(!soma_sparse->has_metadata("md"));
        REQUIRE(soma_sparse->metadata_num() == 2);
    }
}
