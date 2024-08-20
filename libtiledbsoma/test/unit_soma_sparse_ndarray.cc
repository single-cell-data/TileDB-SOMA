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
#define DIM_MAX 1000

TEST_CASE("SOMASparseNDArray: basic") {
    int64_t dim_max = 1000;
    bool use_current_domains[] = {false, true};
    for (bool use_current_domain : use_current_domains) {
        auto ctx = std::make_shared<SOMAContext>();
        std::string uri = "mem://unit-test-sparse-ndarray-basic";

        REQUIRE(!SOMASparseNDArray::exists(uri, ctx));

        auto index_columns = helper::create_column_index_info(
            dim_max, use_current_domain);
        SOMASparseNDArray::create(
            uri,
            "l",
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
        REQUIRE(soma_sparse->soma_data_type() == "l");
        auto schema = soma_sparse->tiledb_schema();
        REQUIRE(schema->has_attribute("soma_data"));
        REQUIRE(schema->array_type() == TILEDB_SPARSE);
        REQUIRE(schema->domain().has_dimension("soma_dim_0"));
        REQUIRE(soma_sparse->ndim() == 1);
        REQUIRE(soma_sparse->nnz() == 0);

        if (use_current_domain) {
            REQUIRE(soma_sparse->shape() == std::vector<int64_t>{dim_max + 1});
        } else {
            REQUIRE(
                soma_sparse->maxshape() == std::vector<int64_t>{dim_max + 1});
        }

        soma_sparse->close();

        std::vector<int64_t> d0(10);
        for (int j = 0; j < 10; j++)
            d0[j] = j;
        std::vector<int> a0(10, 1);

        soma_sparse->open(OpenMode::write);
        soma_sparse->set_column_data("soma_data", a0.size(), a0.data());
        soma_sparse->set_column_data("soma_dim_0", d0.size(), d0.data());
        soma_sparse->write();
        soma_sparse->close();

        soma_sparse->open(OpenMode::read);
        while (auto batch = soma_sparse->read_next()) {
            auto arrbuf = batch.value();
            auto d0span = arrbuf->at("soma_dim_0")->data<int64_t>();
            auto a0span = arrbuf->at("soma_data")->data<int>();
            REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
            REQUIRE(a0 == std::vector<int>(a0span.begin(), a0span.end()));
        }

        // TODO on a subsequent PR if use_current_domain:
        // * test write out of bounds, including assertion of exception type
        // * test resize
        // * test write within new bounds

        soma_sparse->close();
    }
}

TEST_CASE("SOMASparseNDArray: platform_config") {
    int64_t dim_max = 1000;
    bool use_current_domains[] = {false, true};
    for (bool use_current_domain : use_current_domains) {
        auto ctx = std::make_shared<SOMAContext>();
        std::string uri = "mem://unit-test-dataframe-platform-config";

        PlatformConfig platform_config;
        platform_config.sparse_nd_array_dim_zstd_level = 6;

        auto index_columns = helper::create_column_index_info(
            dim_max, use_current_domain);
        SOMASparseNDArray::create(
            uri,
            "l",
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            platform_config);

        auto soma_dataframe = SOMASparseNDArray::open(uri, OpenMode::read, ctx);
        auto dim_filter = soma_dataframe->tiledb_schema()
                              ->domain()
                              .dimension("soma_dim_0")
                              .filter_list()
                              .filter(0);
        REQUIRE(dim_filter.filter_type() == TILEDB_FILTER_ZSTD);
        REQUIRE(dim_filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL) == 6);

        soma_dataframe->close();
    }
}

TEST_CASE("SOMASparseNDArray: metadata") {
    int64_t dim_max = 1000;
    bool use_current_domains[] = {false, true};
    for (bool use_current_domain : use_current_domains) {
        auto ctx = std::make_shared<SOMAContext>();

        std::string uri = "mem://unit-test-sparse-ndarray";

        auto index_columns = helper::create_column_index_info(
            dim_max, use_current_domain);
        SOMASparseNDArray::create(
            uri,
            "l",
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
