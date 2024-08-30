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

// This is a keystroke-reduction fixture for some similar unit-test cases For
// convenience there are dims/attrs of type int64, uint32, and string. (Feel
// free to add more types.) The main value-adds of this fixture are (a) simple
// keystroke-reduction; (b) you get to pick which ones are the dim(s) and which
// are the attr(s).
struct VariouslyIndexedDataFrameFixture {
    std::shared_ptr<SOMAContext> ctx_;
    std::string uri_;
    bool use_current_domain_;

    //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Using Catch2's TEST_CASE_METHOD we can't pass constructor args.
    // This is a call-after-construction method.
    void set_up(std::shared_ptr<SOMAContext> ctx, std::string uri) {
        ctx_ = ctx;
        uri_ = uri;
    }

    //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Helpers for setting up dim/attr configs and data
    static const inline int64_t i64_dim_max = 100;
    static const inline int64_t u32_dim_max = 10000;
    static const inline int64_t str_dim_max = 0;  // not used for string dims

    static const inline std::string i64_name = "soma_joinid";
    static const inline std::string u32_name = "myuint32";
    static const inline std::string str_name = "mystring";

    tiledb_datatype_t i64_datatype = TILEDB_INT64;
    tiledb_datatype_t u32_datatype = TILEDB_UINT32;
    tiledb_datatype_t str_datatype = TILEDB_STRING_ASCII;

    std::string i64_arrow_format = helper::to_arrow_format(i64_datatype);
    std::string u32_arrow_format = helper::to_arrow_format(u32_datatype);
    std::string attr_1_arrow_format = helper::to_arrow_format(str_datatype);

    helper::DimInfo i64_dim_info(bool use_current_domain) {
        return helper::DimInfo(
            {.name = i64_name,
             .tiledb_datatype = i64_datatype,
             .dim_max = i64_dim_max,
             .use_current_domain = use_current_domain});
    }
    helper::DimInfo u32_dim_info(bool use_current_domain) {
        return helper::DimInfo(
            {.name = u32_name,
             .tiledb_datatype = u32_datatype,
             .dim_max = u32_dim_max,
             .use_current_domain = use_current_domain});
    }
    helper::DimInfo str_dim_info(bool use_current_domain) {
        return helper::DimInfo(
            {.name = str_name,
             .tiledb_datatype = str_datatype,
             .dim_max = str_dim_max,
             .use_current_domain = use_current_domain});
    }

    helper::AttrInfo i64_attr_info(std::string name = i64_name) {
        return helper::AttrInfo(
            {.name = name, .tiledb_datatype = i64_datatype});
    }
    helper::AttrInfo u32_attr_info() {
        return helper::AttrInfo(
            {.name = u32_name, .tiledb_datatype = u32_datatype});
    }
    helper::AttrInfo str_attr_info() {
        return helper::AttrInfo(
            {.name = str_name, .tiledb_datatype = str_datatype});
    }

    std::vector<int64_t> make_i64_data() {
        return std::vector<int64_t>({1, 2});
    }
    std::vector<uint32_t> make_u32_data() {
        return std::vector<uint32_t>({1234, 5678});
    }
    std::vector<std::string> make_str_data() {
        return std::vector<std::string>({"apple", "bat"});
    }

    //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Helper methods for create/open/write/etc.

    void create(
        const std::vector<helper::DimInfo>& dim_infos,
        const std::vector<helper::AttrInfo>& attr_infos) {
        auto [schema, index_columns] =
            helper::create_arrow_schema_and_index_columns(
                dim_infos, attr_infos);
        SOMADataFrame::create(
            uri_,
            std::move(schema),
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx_);
    }

    void create(
        const std::vector<helper::DimInfo>& dim_infos,
        const std::vector<helper::AttrInfo>& attr_infos,
        const PlatformConfig& platform_config,
        std::optional<TimestampRange> timestamp_range = std::nullopt) {
        auto [schema, index_columns] =
            helper::create_arrow_schema_and_index_columns(
                dim_infos, attr_infos);
        SOMADataFrame::create(
            uri_,
            std::move(schema),
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx_,
            platform_config,
            timestamp_range);
    }

    std::unique_ptr<SOMADataFrame> open(
        OpenMode mode,
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<TimestampRange> timestamp_range = std::nullopt) {
        return SOMADataFrame::open(
            uri_,
            mode,
            ctx_,
            {},  // column_names
            result_order,
            timestamp_range);
    }

    void write_generic_data() {
        // No arguments -- for now.
        // In a subsequent PR we'll vary writing in-bounds vs out of bounds.
        auto soma_dataframe = SOMADataFrame::open(uri_, OpenMode::write, ctx_);

        auto i64_data = make_i64_data();
        auto u32_data = make_u32_data();
        auto str_data = make_str_data();

        soma_dataframe->set_column_data(
            i64_name, i64_data.size(), i64_data.data());
        soma_dataframe->set_column_data(
            str_name, str_data.size(), str_data.data());
        soma_dataframe->set_column_data(
            u32_name, u32_data.size(), u32_data.data());
        soma_dataframe->write();

        soma_dataframe->close();
    }
};

TEST_CASE_METHOD(VariouslyIndexedDataFrameFixture, "SOMADataFrame: basic") {
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        set_up(
            std::make_shared<SOMAContext>(), "mem://unit-test-dataframe-basic");

        std::vector<helper::DimInfo> dim_infos(
            {i64_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

        REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

        create(dim_infos, attr_infos);

        REQUIRE(SOMADataFrame::exists(uri_, ctx_));
        REQUIRE(!SOMASparseNDArray::exists(uri_, ctx_));
        REQUIRE(!SOMADenseNDArray::exists(uri_, ctx_));

        auto soma_dataframe = open(OpenMode::read);
        REQUIRE(soma_dataframe->uri() == uri_);
        REQUIRE(soma_dataframe->ctx() == ctx_);
        REQUIRE(soma_dataframe->type() == "SOMADataFrame");
        std::vector<std::string> expected_index_column_names = {
            dim_infos[0].name};
        REQUIRE(
            soma_dataframe->index_column_names() ==
            expected_index_column_names);
        REQUIRE(soma_dataframe->count() == 0);
        soma_dataframe->close();

        std::vector<int64_t> d0(10);
        for (int j = 0; j < 10; j++)
            d0[j] = j;
        std::vector<int> a0(10, 1);

        soma_dataframe = open(OpenMode::write);
        soma_dataframe->set_column_data(
            dim_infos[0].name, d0.size(), d0.data());
        soma_dataframe->set_column_data(
            attr_infos[0].name, a0.size(), a0.data());
        soma_dataframe->write();
        soma_dataframe->close();

        soma_dataframe = open(OpenMode::read);
        while (auto batch = soma_dataframe->read_next()) {
            auto arrbuf = batch.value();
            auto d0span = arrbuf->at(dim_infos[0].name)->data<int64_t>();
            auto a0span = arrbuf->at(attr_infos[0].name)->data<int>();
            REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
            REQUIRE(a0 == std::vector<int>(a0span.begin(), a0span.end()));
        }
        soma_dataframe->close();

        auto soma_object = SOMAObject::open(uri_, OpenMode::read, ctx_);
        REQUIRE(soma_object->uri() == uri_);
        REQUIRE(soma_object->type() == "SOMADataFrame");
        soma_object->close();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture, "SOMADataFrame: platform_config") {
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

    // TODO this used to be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- filter=" << filter.first;

    SECTION(section.str()) {
        auto use_current_domain = GENERATE(false, true);
        std::ostringstream section2;
        section2 << "- use_current_domain=" << use_current_domain;
        SECTION(section2.str()) {
            set_up(
                std::make_shared<SOMAContext>(),
                "mem://unit-test-dataframe-platform-config");

            PlatformConfig platform_config;
            platform_config.dataframe_dim_zstd_level = 6;
            platform_config.offsets_filters = R"([)" + filter.first + R"(])";
            platform_config.validity_filters = R"([)" + filter.first + R"(])";
            if (filter.second != TILEDB_FILTER_WEBP) {
                platform_config.attrs = R"({"a0": {"filters":[)" +
                                        filter.first + R"(]}})";
            }

            std::vector<helper::DimInfo> dim_infos(
                {i64_dim_info(use_current_domain)});
            std::vector<helper::AttrInfo> attr_infos({i64_attr_info("a0")});

            REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

            create(dim_infos, attr_infos, platform_config);

            auto soma_dataframe = open(OpenMode::read);
            auto sch = soma_dataframe->tiledb_schema();
            REQUIRE(
                sch->offsets_filter_list().filter(0).filter_type() ==
                filter.second);

            REQUIRE(
                sch->validity_filter_list().filter(0).filter_type() ==
                filter.second);

            auto dim_filter = sch->domain()
                                  .dimension(dim_infos[0].name)
                                  .filter_list()
                                  .filter(0);
            REQUIRE(dim_filter.filter_type() == TILEDB_FILTER_ZSTD);
            REQUIRE(
                dim_filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL) == 6);

            if (filter.second != TILEDB_FILTER_WEBP) {
                REQUIRE(
                    sch->attribute(attr_infos[0].name)
                        .filter_list()
                        .filter(0)
                        .filter_type() == filter.second);
            }
            soma_dataframe->close();
        }
    }
}

TEST_CASE_METHOD(VariouslyIndexedDataFrameFixture, "SOMADataFrame: metadata") {
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        set_up(std::make_shared<SOMAContext>(), "mem://unit-test-collection");

        std::vector<helper::DimInfo> dim_infos(
            {i64_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

        REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

        create(dim_infos, attr_infos, PlatformConfig(), TimestampRange(0, 2));

        auto soma_dataframe = open(
            OpenMode::write, ResultOrder::automatic, TimestampRange(1, 1));

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

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture, "SOMADataFrame: bounds-checking") {
    bool use_current_domain = true;
    int old_max = 100;
    int new_max = 200;

    set_up(std::make_shared<SOMAContext>(), "mem://unit-test-bounds-checking");

    std::vector<helper::DimInfo> dim_infos({i64_dim_info(use_current_domain)});
    std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

    REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

    create(dim_infos, attr_infos);

    auto soma_dataframe = open(OpenMode::write);

    std::vector<int64_t> d0({old_max + 1, old_max + 2});
    std::vector<double> a0({1.5, 2.5});
    soma_dataframe->set_column_data(dim_infos[0].name, d0.size(), d0.data());
    soma_dataframe->set_column_data(attr_infos[0].name, a0.size(), a0.data());
    // Writing outside the current domain should fail
    REQUIRE_THROWS(soma_dataframe->write());
    soma_dataframe->close();

    soma_dataframe = open(OpenMode::write);

    soma_dataframe->resize(std::vector<int64_t>({new_max}));
    soma_dataframe->close();

    soma_dataframe = open(OpenMode::write);
    soma_dataframe->set_column_data(dim_infos[0].name, d0.size(), d0.data());
    soma_dataframe->set_column_data(attr_infos[0].name, a0.size(), a0.data());
    // Writing after resize should succeed
    soma_dataframe->write();

    soma_dataframe->close();
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe 1") {
    // LOG_SET_LEVEL("debug");
    auto use_current_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-variant-indexed-dataframe-1");

        std::vector<helper::DimInfo> dim_infos(
            {i64_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos(
            {str_attr_info(), u32_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto soma_dataframe = open(OpenMode::read);

        CurrentDomain current_domain = soma_dataframe->get_current_domain();
        if (!use_current_domain) {
            REQUIRE(current_domain.is_empty());
        } else {
            REQUIRE(!current_domain.is_empty());
            REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
            NDRectangle ndrect = current_domain.ndrectangle();

            std::array<int64_t, 2> i64_range = ndrect.range<int64_t>(
                dim_infos[0].name);
            REQUIRE(i64_range[0] == (int64_t)0);
            REQUIRE(i64_range[1] == (int64_t)dim_infos[0].dim_max);
        }

        soma_dataframe->close();

        write_generic_data();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe 2") {
    // LOG_SET_LEVEL("debug");
    auto use_current_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-variant-indexed-dataframe-2");

        std::vector<helper::DimInfo> dim_infos(
            {i64_dim_info(use_current_domain),
             u32_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos({str_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto soma_dataframe = open(OpenMode::read);

        CurrentDomain current_domain = soma_dataframe->get_current_domain();
        if (!use_current_domain) {
            REQUIRE(current_domain.is_empty());
        } else {
            REQUIRE(!current_domain.is_empty());
            REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
            NDRectangle ndrect = current_domain.ndrectangle();

            std::array<int64_t, 2> i64_range = ndrect.range<int64_t>(
                dim_infos[0].name);
            REQUIRE(i64_range[0] == (int64_t)0);
            REQUIRE(i64_range[1] == (int64_t)dim_infos[0].dim_max);

            std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(
                dim_infos[0].name);
            REQUIRE(u32_range[0] == (uint32_t)0);
            REQUIRE(u32_range[1] == (uint32_t)dim_infos[0].dim_max);
        }

        soma_dataframe->close();

        // Write
        write_generic_data();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe 3") {
    // LOG_SET_LEVEL("debug");
    auto use_current_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-variant-indexed-dataframe-3");

        std::vector<helper::DimInfo> dim_infos(
            {i64_dim_info(use_current_domain),
             str_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto soma_dataframe = open(OpenMode::read);

        CurrentDomain current_domain = soma_dataframe->get_current_domain();
        if (!use_current_domain) {
            REQUIRE(current_domain.is_empty());
        } else {
            REQUIRE(!current_domain.is_empty());
            REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
            NDRectangle ndrect = current_domain.ndrectangle();

            std::array<int64_t, 2> i64_range = ndrect.range<int64_t>(
                dim_infos[0].name);
            REQUIRE(i64_range[0] == (int64_t)0);
            REQUIRE(i64_range[1] == (int64_t)dim_infos[0].dim_max);

            std::array<std::string, 2> str_range = ndrect.range<std::string>(
                dim_infos[1].name);
            // Can we write empty strings in this range?
            REQUIRE(str_range[0] <= "");
            REQUIRE(str_range[1] >= "");
            // Can we write ASCII values in this range?
            REQUIRE(str_range[0] < " ");
            REQUIRE(str_range[1] > "~");
        }

        soma_dataframe->close();

        // Write
        write_generic_data();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe 4") {
    // LOG_SET_LEVEL("debug");
    auto use_current_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-variant-indexed-dataframe-4");

        std::vector<helper::DimInfo> dim_infos(
            {str_dim_info(use_current_domain),
             u32_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos({i64_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto soma_dataframe = open(OpenMode::read);

        CurrentDomain current_domain = soma_dataframe->get_current_domain();
        if (!use_current_domain) {
            REQUIRE(current_domain.is_empty());
        } else {
            REQUIRE(!current_domain.is_empty());
            REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
            NDRectangle ndrect = current_domain.ndrectangle();

            std::array<std::string, 2> str_range = ndrect.range<std::string>(
                dim_infos[0].name);
            // Can we write empty strings in this range?
            REQUIRE(str_range[0] <= "");
            REQUIRE(str_range[1] >= "");
            // Can we write ASCII values in this range?
            REQUIRE(str_range[0] < " ");
            REQUIRE(str_range[1] > "~");

            std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(
                dim_infos[1].name);
            REQUIRE(u32_range[0] == (uint32_t)0);
            REQUIRE(u32_range[1] == (uint32_t)dim_infos[1].dim_max);
        }

        soma_dataframe->close();

        // Write
        write_generic_data();
    }
}
