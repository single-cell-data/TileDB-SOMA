/**
 * @file   unit_soma_dataframe.cc
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
 * This file manages unit tests for the SOMADataFrame class
 */

#include "common.h"

#define DIM_MAX 1000

TEST_CASE("SOMADataFrame: basic") {
    int64_t dim_max = 1000;
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();
        std::string uri = "mem://unit-test-dataframe-basic";

        REQUIRE(!SOMADataFrame::exists(uri, ctx));

        auto [schema, index_columns] =
            helper::create_arrow_schema_and_index_columns(
                dim_max, use_current_domain);
        SOMADataFrame::create(
            uri,
            std::move(schema),
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx);

        REQUIRE(SOMADataFrame::exists(uri, ctx));
        REQUIRE(!SOMASparseNDArray::exists(uri, ctx));
        REQUIRE(!SOMADenseNDArray::exists(uri, ctx));

        auto soma_dataframe = SOMADataFrame::open(uri, OpenMode::read, ctx);
        REQUIRE(soma_dataframe->uri() == uri);
        REQUIRE(soma_dataframe->ctx() == ctx);
        REQUIRE(soma_dataframe->type() == "SOMADataFrame");
        std::vector<std::string> expected_index_column_names = {"d0"};
        REQUIRE(
            soma_dataframe->index_column_names() ==
            expected_index_column_names);
        REQUIRE(soma_dataframe->count() == 0);
        soma_dataframe->close();

        std::vector<int64_t> d0(10);
        for (int j = 0; j < 10; j++)
            d0[j] = j;
        std::vector<int> a0(10, 1);

        soma_dataframe = SOMADataFrame::open(uri, OpenMode::write, ctx);
        soma_dataframe->set_column_data("a0", a0.size(), a0.data());
        soma_dataframe->set_column_data("d0", d0.size(), d0.data());
        soma_dataframe->write();
        soma_dataframe->close();

        soma_dataframe = SOMADataFrame::open(uri, OpenMode::read, ctx);
        while (auto batch = soma_dataframe->read_next()) {
            auto arrbuf = batch.value();
            auto d0span = arrbuf->at("d0")->data<int64_t>();
            auto a0span = arrbuf->at("a0")->data<int>();
            REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
            REQUIRE(a0 == std::vector<int>(a0span.begin(), a0span.end()));
        }
        soma_dataframe->close();

        auto soma_object = SOMAObject::open(uri, OpenMode::read, ctx);
        REQUIRE(soma_object->uri() == uri);
        REQUIRE(soma_object->type() == "SOMADataFrame");
        soma_object->close();
    }
}

TEST_CASE("SOMADataFrame: platform_config") {
    std::pair<std::string, tiledb_filter_type_t> filter = GENERATE(
        std::make_pair(
            R"({"name": "GZIP", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_GZIP),
        std::make_pair(
            R"({"name": "ZSTD", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_ZSTD),
        std::make_pair(
            R"({"name": "LZ4", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_LZ4),
        std::make_pair(
            R"({"name": "BZIP2", "COMPRESSION_LEVEL": 3})",
            TILEDB_FILTER_BZIP2),
        std::make_pair(
            R"({"name": "RLE", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_RLE),
        std::make_pair(
            R"({"name": "DICTIONARY_ENCODING", "COMPRESSION_LEVEL": 3})",
            TILEDB_FILTER_DICTIONARY),
        std::make_pair(
            R"({"name": "BIT_WIDTH_REDUCTION", "BIT_WIDTH_MAX_WINDOW": 3})",
            TILEDB_FILTER_BIT_WIDTH_REDUCTION),
        std::make_pair(
            R"({"name": "POSITIVE_DELTA", "POSITIVE_DELTA_MAX_WINDOW": 3})",
            TILEDB_FILTER_POSITIVE_DELTA),
        std::make_pair(
            R"({"name": "DELTA", "COMPRESSION_LEVEL": 3, "COMPRESSION_REINTERPRET_DATATYPE": 1})",
            TILEDB_FILTER_DELTA),
        std::make_pair(
            R"({"name": "DOUBLE_DELTA", "COMPRESSION_LEVEL": 3, "COMPRESSION_REINTERPRET_DATATYPE": 1})",
            TILEDB_FILTER_DOUBLE_DELTA),
        std::make_pair(
            R"({"name": "SCALE_FLOAT", "SCALE_FLOAT_FACTOR": 1, "SCALE_FLOAT_OFFSET": 0, "SCALE_FLOAT_BYTEWIDTH": 8})",
            TILEDB_FILTER_SCALE_FLOAT),
        std::make_pair(
            R"({"name": "WEBP", "WEBP_INPUT_FORMAT": 1, "WEBP_QUALITY": 25.5,
            "WEBP_LOSSLESS": 0})",
            TILEDB_FILTER_WEBP),
        std::make_pair(R"("CHECKSUM_MD5")", TILEDB_FILTER_CHECKSUM_MD5),
        std::make_pair(R"("CHECKSUM_SHA256")", TILEDB_FILTER_CHECKSUM_SHA256),
        std::make_pair(R"("XOR")", TILEDB_FILTER_XOR),
        std::make_pair(R"("BITSHUFFLE")", TILEDB_FILTER_BITSHUFFLE),
        std::make_pair(R"("BYTESHUFFLE")", TILEDB_FILTER_BYTESHUFFLE),
        std::make_pair(R"("NOOP")", TILEDB_FILTER_NONE));

    // TODO this use to be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- filter=" << filter.first;

    SECTION(section.str()) {
        int64_t dim_max = 1000;
        auto use_current_domain = GENERATE(false, true);
        // TODO this could be formatted with fmt::format which is part of
        // internal header spd/log/fmt/fmt.h and should not be used. In C++20,
        // this can be replaced with std::format.
        std::ostringstream section2;
        section2 << "- use_current_domain=" << use_current_domain;
        SECTION(section2.str()) {
            auto ctx = std::make_shared<SOMAContext>();
            std::string uri = "mem://unit-test-dataframe-platform-config";

            PlatformConfig platform_config;
            platform_config.dataframe_dim_zstd_level = 6;
            platform_config.offsets_filters = R"([)" + filter.first + R"(])";
            platform_config.validity_filters = R"([)" + filter.first + R"(])";
            if (filter.second != TILEDB_FILTER_WEBP) {
                platform_config.attrs = R"({"a0": {"filters":[)" +
                                        filter.first + R"(]}})";
            }

            auto [schema, index_columns] =
                helper::create_arrow_schema_and_index_columns(
                    dim_max, use_current_domain);
            SOMADataFrame::create(
                uri,
                std::move(schema),
                ArrowTable(
                    std::move(index_columns.first),
                    std::move(index_columns.second)),
                ctx,
                platform_config);

            auto soma_dataframe = SOMADataFrame::open(uri, OpenMode::read, ctx);
            auto sch = soma_dataframe->tiledb_schema();
            REQUIRE(
                sch->offsets_filter_list().filter(0).filter_type() ==
                filter.second);

            REQUIRE(
                sch->validity_filter_list().filter(0).filter_type() ==
                filter.second);

            auto dim_filter = sch->domain()
                                  .dimension("d0")
                                  .filter_list()
                                  .filter(0);
            REQUIRE(dim_filter.filter_type() == TILEDB_FILTER_ZSTD);
            REQUIRE(
                dim_filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL) == 6);

            if (filter.second != TILEDB_FILTER_WEBP) {
                REQUIRE(
                    sch->attribute("a0")
                        .filter_list()
                        .filter(0)
                        .filter_type() == filter.second);
            }
            soma_dataframe->close();
        }
    }
}

TEST_CASE("SOMADataFrame: metadata") {
    int64_t dim_max = 1000;
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();
        std::string uri = "mem://unit-test-collection";

        auto [schema, index_columns] =
            helper::create_arrow_schema_and_index_columns(
                dim_max, use_current_domain);
        SOMADataFrame::create(
            uri,
            std::move(schema),
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            PlatformConfig(),
            TimestampRange(0, 2));

        auto soma_dataframe = SOMADataFrame::open(
            uri,
            OpenMode::write,
            ctx,
            {},
            ResultOrder::automatic,
            TimestampRange(1, 1));

        int32_t val = 100;
        soma_dataframe->set_metadata("md", TILEDB_INT32, 1, &val);
        soma_dataframe->close();

        // Read metadata
        soma_dataframe->open(OpenMode::read, TimestampRange(0, 2));
        REQUIRE(soma_dataframe->metadata_num() == 3);
        REQUIRE(soma_dataframe->has_metadata("soma_object_type"));
        REQUIRE(soma_dataframe->has_metadata("soma_encoding_version"));
        REQUIRE(soma_dataframe->has_metadata("md"));
        auto mdval = soma_dataframe->get_metadata("md");
        REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
        REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
        soma_dataframe->close();

        // md should not be available at (2, 2)
        soma_dataframe->open(OpenMode::read, TimestampRange(2, 2));
        REQUIRE(soma_dataframe->metadata_num() == 2);
        REQUIRE(soma_dataframe->has_metadata("soma_object_type"));
        REQUIRE(soma_dataframe->has_metadata("soma_encoding_version"));
        REQUIRE(!soma_dataframe->has_metadata("md"));
        soma_dataframe->close();

        // Metadata should also be retrievable in write mode
        soma_dataframe->open(OpenMode::write, TimestampRange(0, 2));
        REQUIRE(soma_dataframe->metadata_num() == 3);
        REQUIRE(soma_dataframe->has_metadata("soma_object_type"));
        REQUIRE(soma_dataframe->has_metadata("soma_encoding_version"));
        REQUIRE(soma_dataframe->has_metadata("md"));
        mdval = soma_dataframe->get_metadata("md");
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

        // Delete and have it reflected when reading metadata while in write
        // mode
        soma_dataframe->delete_metadata("md");
        mdval = soma_dataframe->get_metadata("md");
        REQUIRE(!mdval.has_value());
        soma_dataframe->close();

        // Confirm delete in read mode
        soma_dataframe->open(OpenMode::read, TimestampRange(0, 2));
        REQUIRE(!soma_dataframe->has_metadata("md"));
        REQUIRE(soma_dataframe->metadata_num() == 2);
    }
}

TEST_CASE("SOMADataFrame: bounds-checking") {
    bool use_current_domain = true;
    int old_max = 100;
    int new_max = 200;

    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-bounds-checking";

    auto [schema, index_columns] =
        helper::create_arrow_schema_and_index_columns(100, use_current_domain);

    SOMADataFrame::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx);

    auto soma_dataframe = SOMADataFrame::open(uri, OpenMode::write, ctx);

    std::vector<int64_t> d0({old_max + 1, old_max + 2});
    std::vector<double> a0({1.5, 2.5});
    soma_dataframe->set_column_data("d0", d0.size(), d0.data());
    soma_dataframe->set_column_data("a0", a0.size(), a0.data());
    // Writing outside the current domain should fail
    REQUIRE_THROWS(soma_dataframe->write());
    soma_dataframe->close();

    soma_dataframe = SOMADataFrame::open(uri, OpenMode::write, ctx);
    soma_dataframe->resize(std::vector<int64_t>({new_max}));
    soma_dataframe->close();

    soma_dataframe = SOMADataFrame::open(uri, OpenMode::write, ctx);
    soma_dataframe->set_column_data("d0", d0.size(), d0.data());
    soma_dataframe->set_column_data("a0", a0.size(), a0.data());
    // Writing after resize should succeed
    soma_dataframe->write();

    soma_dataframe->close();
}
