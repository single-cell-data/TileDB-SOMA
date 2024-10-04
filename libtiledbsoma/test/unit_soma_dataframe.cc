/**
 * @file   unit_sdf.cc
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

const int64_t SOMA_JOINID_DIM_MAX = 99;
const int64_t SOMA_JOINID_RESIZE_DIM_MAX = 199;

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
    static const inline int64_t i64_dim_max = SOMA_JOINID_DIM_MAX;
    static const inline int64_t u32_dim_max = 9999;
    static const inline int64_t str_dim_max = 0;  // not used for string dims

    static const inline std::string i64_name = "soma_joinid";
    static const inline std::string u32_name = "myuint32";
    static const inline std::string str_name = "mystring";

    tiledb_datatype_t i64_datatype = TILEDB_INT64;
    tiledb_datatype_t u32_datatype = TILEDB_UINT32;
    tiledb_datatype_t str_datatype = TILEDB_STRING_ASCII;

    std::string i64_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        i64_datatype);
    std::string u32_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        u32_datatype);
    std::string attr_1_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        str_datatype);

    helper::DimInfo i64_dim_info(bool use_current_domain) {
        return helper::DimInfo(
            {.name = i64_name,
             .tiledb_datatype = i64_datatype,
             .dim_max = i64_dim_max,
             .string_lo = "N/A",
             .string_hi = "N/A",
             .use_current_domain = use_current_domain});
    }
    helper::DimInfo u32_dim_info(bool use_current_domain) {
        return helper::DimInfo(
            {.name = u32_name,
             .tiledb_datatype = u32_datatype,
             .dim_max = u32_dim_max,
             .string_lo = "N/A",
             .string_hi = "N/A",
             .use_current_domain = use_current_domain});
    }
    helper::DimInfo str_dim_info(
        bool use_current_domain, std::string string_lo, std::string string_hi) {
        return helper::DimInfo(
            {.name = str_name,
             .tiledb_datatype = str_datatype,
             .dim_max = str_dim_max,
             .string_lo = string_lo,
             .string_hi = string_hi,
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

    void write_sjid_u32_str_data_from(int64_t sjid_base) {
        auto sdf = SOMADataFrame::open(uri_, OpenMode::write, ctx_);

        auto i64_data = std::vector<int64_t>({sjid_base + 1, sjid_base + 2});

        auto u32_data = std::vector<uint32_t>({1234, 5678});

        // We like to think we're writing an array of strings ...
        auto strings = std::vector<std::string>({"apple", "bat"});
        // ... but really we're writing an array of characters along
        // with offsets data.
        //
        // It would be possible here to just hard-code a string "applebat" and
        // an offsets array {0, 5, 8}. The following bits simply automate that.
        std::string char_data("");
        std::vector<uint64_t> char_offsets(0);
        uint64_t offset = 0;
        for (auto e : strings) {
            char_data += e;
            char_offsets.push_back(offset);
            offset += e.size();
        }
        char_offsets.push_back(offset);

        sdf->set_column_data(i64_name, i64_data.size(), i64_data.data());
        sdf->set_column_data(
            str_name, strings.size(), char_data.data(), char_offsets.data());
        sdf->set_column_data(u32_name, u32_data.size(), u32_data.data());
        sdf->write();

        sdf->close();
    }
};

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: basic",
    "[SOMADataFrame]") {
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

        auto sdf = open(OpenMode::read);
        REQUIRE(sdf->uri() == uri_);
        REQUIRE(sdf->ctx() == ctx_);
        REQUIRE(sdf->type() == "SOMADataFrame");
        std::vector<std::string> expected_index_column_names = {
            dim_infos[0].name};
        REQUIRE(sdf->index_column_names() == expected_index_column_names);
        REQUIRE(sdf->nnz() == 0);
        sdf->close();

        std::vector<int64_t> d0(10);
        for (int j = 0; j < 10; j++)
            d0[j] = j;
        std::vector<uint32_t> a0(10, 1);

        sdf = open(OpenMode::write);
        sdf->set_column_data(dim_infos[0].name, d0.size(), d0.data());
        sdf->set_column_data(attr_infos[0].name, a0.size(), a0.data());
        sdf->write();
        sdf->close();

        sdf = open(OpenMode::read);
        while (auto batch = sdf->read_next()) {
            auto arrbuf = batch.value();
            auto d0span = arrbuf->at(dim_infos[0].name)->data<int64_t>();
            auto a0span = arrbuf->at(attr_infos[0].name)->data<uint32_t>();
            REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
            REQUIRE(a0 == std::vector<uint32_t>(a0span.begin(), a0span.end()));
        }
        sdf->close();

        auto soma_object = SOMAObject::open(uri_, OpenMode::read, ctx_);
        REQUIRE(soma_object->uri() == uri_);
        REQUIRE(soma_object->type() == "SOMADataFrame");
        soma_object->close();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: platform_config",
    "[SOMADataFrame]") {
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

            auto sdf = open(OpenMode::read);
            auto sch = sdf->tiledb_schema();
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
            sdf->close();
        }
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: metadata",
    "[SOMADataFrame]") {
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

        auto sdf = open(
            OpenMode::write, ResultOrder::automatic, TimestampRange(1, 1));

        int32_t val = 100;
        sdf->set_metadata("md", TILEDB_INT32, 1, &val);
        sdf->close();

        // Read metadata
        sdf->open(OpenMode::read, TimestampRange(0, 2));
        REQUIRE(sdf->metadata_num() == 3);
        REQUIRE(sdf->has_metadata("soma_object_type"));
        REQUIRE(sdf->has_metadata("soma_encoding_version"));
        REQUIRE(sdf->has_metadata("md"));
        auto mdval = sdf->get_metadata("md");
        REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
        REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
        sdf->close();

        // md should not be available at (2, 2)
        sdf->open(OpenMode::read, TimestampRange(2, 2));
        REQUIRE(sdf->metadata_num() == 2);
        REQUIRE(sdf->has_metadata("soma_object_type"));
        REQUIRE(sdf->has_metadata("soma_encoding_version"));
        REQUIRE(!sdf->has_metadata("md"));
        sdf->close();

        // Metadata should also be retrievable in write mode
        sdf->open(OpenMode::write, TimestampRange(0, 2));
        REQUIRE(sdf->metadata_num() == 3);
        REQUIRE(sdf->has_metadata("soma_object_type"));
        REQUIRE(sdf->has_metadata("soma_encoding_version"));
        REQUIRE(sdf->has_metadata("md"));
        mdval = sdf->get_metadata("md");
        REQUIRE(
            *((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

        // Delete and have it reflected when reading metadata while in write
        // mode
        sdf->delete_metadata("md");
        mdval = sdf->get_metadata("md");
        REQUIRE(!mdval.has_value());
        sdf->close();

        // Confirm delete in read mode
        sdf->open(OpenMode::read, TimestampRange(0, 2));
        REQUIRE(!sdf->has_metadata("md"));
        REQUIRE(sdf->metadata_num() == 2);
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: bounds-checking",
    "[SOMADataFrame]") {
    bool use_current_domain = true;
    int old_max = SOMA_JOINID_DIM_MAX;
    int new_max = SOMA_JOINID_RESIZE_DIM_MAX;

    set_up(std::make_shared<SOMAContext>(), "mem://unit-test-bounds-checking");

    std::vector<helper::DimInfo> dim_infos({i64_dim_info(use_current_domain)});
    std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

    REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

    create(dim_infos, attr_infos);

    auto sdf = open(OpenMode::write);

    std::vector<int64_t> d0({old_max + 1, old_max + 2});
    std::vector<double> a0({1.5, 2.5});
    sdf->set_column_data(dim_infos[0].name, d0.size(), d0.data());
    sdf->set_column_data(attr_infos[0].name, a0.size(), a0.data());
    // Writing outside the current domain should fail
    REQUIRE_THROWS(sdf->write());
    sdf->close();

    sdf = open(OpenMode::write);
    sdf->resize_soma_joinid_shape(int64_t{new_max}, "testing");
    sdf->close();

    sdf = open(OpenMode::write);
    sdf->set_column_data(dim_infos[0].name, d0.size(), d0.data());
    sdf->set_column_data(attr_infos[0].name, a0.size(), a0.data());
    // Writing after resize should succeed
    sdf->write();

    sdf->close();
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: standard-indexed dataframe dim-sjid attr-str-u32",
    "[SOMADataFrame]") {
    auto use_current_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        std::string suffix = use_current_domain ? "true" : "false";
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-variant-indexed-dataframe-1-" + suffix);

        std::vector<helper::DimInfo> dim_infos(
            {i64_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos(
            {str_attr_info(), u32_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::read);

        CurrentDomain current_domain = sdf->get_current_domain_for_test();
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

        // Check shape before resize
        int64_t expect = dim_infos[0].dim_max + 1;
        std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);

        REQUIRE(sdf->nnz() == 0);

        sdf->close();

        // Write data
        write_sjid_u32_str_data_from(0);

        // Check shape after write
        sdf->open(OpenMode::read);

        REQUIRE(sdf->nnz() == 2);

        expect = dim_infos[0].dim_max + 1;
        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);

        // Check domainish accessors before resize
        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<int64_t> ned_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                non_empty_domain, "soma_joinid");

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<int64_t> dom_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                soma_domain, "soma_joinid");

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<int64_t> maxdom_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                soma_maxdomain, "soma_joinid");

        REQUIRE(ned_sjid == std::vector<int64_t>({1, 2}));

        REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));

        REQUIRE(maxdom_sjid.size() == 2);
        REQUIRE(maxdom_sjid[0] == 0);
        if (!use_current_domain) {
            REQUIRE(maxdom_sjid[1] == 99);
        } else {
            REQUIRE(maxdom_sjid[1] > 2000000000);
        }

        sdf->close();

        REQUIRE(sdf->nnz() == 2);
        write_sjid_u32_str_data_from(8);
        REQUIRE(sdf->nnz() == 4);

        // Resize
        auto new_shape = int64_t{SOMA_JOINID_RESIZE_DIM_MAX + 1};

        if (!use_current_domain) {
            // Domain is already set. The domain (not current domain but domain)
            // is immutable. All we can do is check for:
            // * throw on write beyond domain
            // * throw on an attempt to resize.
            REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

            sdf = open(OpenMode::write);
            // Array not resizeable if it has not already been sized
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

        } else {
            // Expect throw on write beyond current domain before resize
            REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

            // Check shape after write
            sdf = open(OpenMode::read);
            expect = dim_infos[0].dim_max + 1;

            std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(actual.has_value());
            REQUIRE(actual.value() == expect);
            sdf->close();

            sdf = open(OpenMode::read);
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

            sdf = open(OpenMode::write);
            sdf->resize_soma_joinid_shape(new_shape, "testing");
            sdf->close();

            // Check shape after resize
            sdf = open(OpenMode::read);
            expect = SOMA_JOINID_RESIZE_DIM_MAX + 1;
            actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(actual.has_value());
            REQUIRE(actual.value() == expect);

            sdf->close();

            // Implicitly we expect no throw
            write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX);
        }

        // Check domainish accessors after resize
        sdf->open(OpenMode::read);

        non_empty_domain = sdf->get_non_empty_domain();
        ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            non_empty_domain, "soma_joinid");

        soma_domain = sdf->get_soma_domain();
        dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            soma_domain, "soma_joinid");

        soma_maxdomain = sdf->get_soma_maxdomain();
        maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<
            int64_t>(soma_maxdomain, "soma_joinid");

        if (!use_current_domain) {
            REQUIRE(ned_sjid == std::vector<int64_t>({1, 10}));
            REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));
            REQUIRE(maxdom_sjid == std::vector<int64_t>({0, 99}));
        } else {
            REQUIRE(ned_sjid == std::vector<int64_t>({1, 101}));
            REQUIRE(dom_sjid == std::vector<int64_t>({0, 199}));
            REQUIRE(maxdom_sjid.size() == 2);
            REQUIRE(maxdom_sjid[0] == 0);
            REQUIRE(maxdom_sjid[1] > 2000000000);
        }

        // Check can_resize_soma_joinid_shape
        std::pair<bool, std::string> check = sdf->can_resize_soma_joinid_shape(
            1, "testing");
        if (!use_current_domain) {
            REQUIRE(check.first == false);
            REQUIRE(
                check.second ==
                "testing: dataframe currently has no domain set: please "
                "upgrade the array.");
        } else {
            // Must fail since this is too small.
            REQUIRE(check.first == false);
            REQUIRE(
                check.second ==
                "testing: new soma_joinid shape 1 < existing shape 199");
            check = sdf->can_resize_soma_joinid_shape(
                SOMA_JOINID_RESIZE_DIM_MAX + 1, "testing");
            REQUIRE(check.first == true);
            REQUIRE(check.second == "");
        }

        sdf->close();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe dim-u32-sjid attr-str",
    "[SOMADataFrame]") {
    auto use_current_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        std::string suffix = use_current_domain ? "true" : "false";
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-variant-indexed-dataframe-2-" + suffix);

        std::vector<helper::DimInfo> dim_infos(
            {u32_dim_info(use_current_domain),
             i64_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos({str_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::read);

        CurrentDomain current_domain = sdf->get_current_domain_for_test();
        if (!use_current_domain) {
            REQUIRE(current_domain.is_empty());
        } else {
            REQUIRE(!current_domain.is_empty());
            REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
            NDRectangle ndrect = current_domain.ndrectangle();

            std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(
                dim_infos[0].name);
            REQUIRE(u32_range[0] == (uint32_t)0);
            REQUIRE(u32_range[1] == (uint32_t)dim_infos[0].dim_max);

            std::array<int64_t, 2> i64_range = ndrect.range<int64_t>(
                dim_infos[1].name);
            REQUIRE(i64_range[0] == (int64_t)0);
            REQUIRE(i64_range[1] == (int64_t)dim_infos[1].dim_max);
        }

        // Check shape before write
        int64_t expect = dim_infos[1].dim_max + 1;
        std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);

        sdf->close();

        REQUIRE(sdf->nnz() == 0);

        // Write
        write_sjid_u32_str_data_from(0);

        REQUIRE(sdf->nnz() == 2);
        write_sjid_u32_str_data_from(8);
        REQUIRE(sdf->nnz() == 4);

        // Check domainish accessors before resize
        sdf->open(OpenMode::read);

        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<int64_t> ned_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                non_empty_domain, "soma_joinid");
        std::vector<uint32_t> ned_u32 =
            ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(
                non_empty_domain, "myuint32");

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<int64_t> dom_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                soma_domain, "soma_joinid");
        std::vector<uint32_t> dom_u32 =
            ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(
                soma_domain, "myuint32");

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<int64_t> maxdom_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                soma_maxdomain, "soma_joinid");
        std::vector<uint32_t> maxdom_u32 =
            ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(
                soma_maxdomain, "myuint32");

        REQUIRE(ned_sjid == std::vector<int64_t>({1, 10}));
        REQUIRE(ned_u32 == std::vector<uint32_t>({1234, 5678}));

        REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));
        REQUIRE(dom_u32 == std::vector<uint32_t>({0, 9999}));

        REQUIRE(maxdom_sjid.size() == 2);
        REQUIRE(maxdom_u32.size() == 2);

        REQUIRE(maxdom_u32[0] == 0);
        if (!use_current_domain) {
            REQUIRE(maxdom_u32[1] == 9999);
        } else {
            REQUIRE(maxdom_u32[1] > 2000000000);
        }

        sdf->close();

        // Check shape after write
        sdf = open(OpenMode::read);
        expect = dim_infos[1].dim_max + 1;
        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);
        sdf->close();

        // Resize
        auto new_shape = int64_t{SOMA_JOINID_RESIZE_DIM_MAX + 1};

        if (!use_current_domain) {
            // Domain is already set. The domain (not current domain but domain)
            // is immutable. All we can do is check for:
            // * throw on write beyond domain
            // * throw on an attempt to resize.
            REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

            sdf = open(OpenMode::write);
            // Array not resizeable if it has not already been sized
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

        } else {
            // Expect throw on write beyond current domain before resize
            REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

            // Check shape after write
            sdf = open(OpenMode::read);
            expect = dim_infos[1].dim_max + 1;
            std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(actual.has_value());
            REQUIRE(actual.value() == expect);
            sdf->close();

            sdf = open(OpenMode::read);
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

            sdf = open(OpenMode::write);
            sdf->resize_soma_joinid_shape(new_shape, "testing");
            sdf->close();

            // Check shape after resize
            sdf = open(OpenMode::read);
            expect = SOMA_JOINID_RESIZE_DIM_MAX + 1;
            actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(actual.has_value());
            REQUIRE(actual.value() == expect);
            sdf->close();

            // Implicitly we expect no throw
            write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX);
        }

        // Check domainish accessors after resize
        sdf->open(OpenMode::read);

        non_empty_domain = sdf->get_non_empty_domain();
        ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            non_empty_domain, "soma_joinid");
        ned_u32 = ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(
            non_empty_domain, "myuint32");

        soma_domain = sdf->get_soma_domain();
        dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            soma_domain, "soma_joinid");
        dom_u32 = ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(
            soma_domain, "myuint32");

        soma_maxdomain = sdf->get_soma_maxdomain();
        maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<
            int64_t>(soma_maxdomain, "soma_joinid");
        maxdom_u32 = ArrowAdapter::get_table_non_string_column_by_name<
            uint32_t>(soma_maxdomain, "myuint32");

        if (!use_current_domain) {
            REQUIRE(ned_sjid == std::vector<int64_t>({1, 10}));
            REQUIRE(ned_u32 == std::vector<uint32_t>({1234, 5678}));

            REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));
            REQUIRE(dom_u32 == std::vector<uint32_t>({0, 9999}));

            REQUIRE(maxdom_sjid == std::vector<int64_t>({0, 99}));
            REQUIRE(maxdom_u32 == std::vector<uint32_t>({0, 9999}));

        } else {
            REQUIRE(ned_sjid == std::vector<int64_t>({1, 101}));
            REQUIRE(ned_u32 == std::vector<uint32_t>({1234, 5678}));

            REQUIRE(dom_sjid == std::vector<int64_t>({0, 199}));
            REQUIRE(dom_u32 == std::vector<uint32_t>({0, 9999}));

            REQUIRE(maxdom_sjid.size() == 2);
            REQUIRE(maxdom_sjid[0] == 0);
            REQUIRE(maxdom_sjid[1] > 2000000000);

            REQUIRE(maxdom_u32.size() == 2);
            REQUIRE(maxdom_u32[0] == 0);
            REQUIRE(maxdom_u32[1] > 2000000000);
        }

        // Check can_resize_soma_joinid_shape
        std::pair<bool, std::string> check = sdf->can_resize_soma_joinid_shape(
            1, "testing");
        if (!use_current_domain) {
            REQUIRE(check.first == false);
            REQUIRE(
                check.second ==
                "testing: dataframe currently has no domain set: please "
                "upgrade the array.");
        } else {
            // Must fail since this is too small.
            REQUIRE(check.first == false);
            REQUIRE(
                check.second ==
                "testing: new soma_joinid shape 1 < existing shape 199");
            check = sdf->can_resize_soma_joinid_shape(
                SOMA_JOINID_RESIZE_DIM_MAX + 1, "testing");
            REQUIRE(check.first == true);
            REQUIRE(check.second == "");
        }

        sdf->close();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe dim-sjid-str attr-u32",
    "[SOMADataFrame]") {
    auto use_current_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        auto specify_domain = GENERATE(false, true);
        std::ostringstream section2;
        section2 << "- specify_domain=" << specify_domain;

        std::string suffix1 = use_current_domain ? "true" : "false";
        std::string suffix2 = specify_domain ? "true" : "false";
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-variant-indexed-dataframe-3-" + suffix1 + "-" +
                suffix2);

        std::string string_lo = specify_domain ? "apple" : "";
        std::string string_hi = specify_domain ? "zebra" : "";
        std::vector<helper::DimInfo> dim_infos(
            {i64_dim_info(use_current_domain),
             str_dim_info(use_current_domain, string_lo, string_hi)});
        std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::read);

        CurrentDomain current_domain = sdf->get_current_domain_for_test();
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
            if (specify_domain) {
                REQUIRE(str_range[0] == dim_infos[1].string_lo);
                REQUIRE(str_range[1] == dim_infos[1].string_hi);
            } else {
                // Can we write empty strings in this range?
                REQUIRE(str_range[0] <= "");
                REQUIRE(str_range[1] >= "");
                // Can we write ASCII values in this range?
                REQUIRE(str_range[0] < " ");
                REQUIRE(str_range[1] > "~");
            }
        }

        // Check shape before write
        int64_t expect = dim_infos[0].dim_max + 1;
        std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);
        sdf->close();

        REQUIRE(sdf->nnz() == 0);

        // Write
        write_sjid_u32_str_data_from(0);

        REQUIRE(sdf->nnz() == 2);
        write_sjid_u32_str_data_from(8);
        REQUIRE(sdf->nnz() == 4);

        // Check shape after write
        sdf = open(OpenMode::read);
        expect = dim_infos[0].dim_max + 1;
        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);

        // Check domainish accessors before resize
        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<int64_t> ned_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                non_empty_domain, "soma_joinid");
        std::vector<std::string>
            ned_str = ArrowAdapter::get_table_string_column_by_name(
                non_empty_domain, "mystring");

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<int64_t> dom_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                soma_domain, "soma_joinid");
        std::vector<std::string>
            dom_str = ArrowAdapter::get_table_string_column_by_name(
                soma_domain, "mystring");

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<int64_t> maxdom_sjid =
            ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                soma_maxdomain, "soma_joinid");
        std::vector<std::string>
            maxdom_str = ArrowAdapter::get_table_string_column_by_name(
                soma_maxdomain, "mystring");

        REQUIRE(ned_sjid == std::vector<int64_t>({1, 10}));
        REQUIRE(ned_str == std::vector<std::string>({"apple", "bat"}));

        REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));

        if (!use_current_domain) {
            REQUIRE(maxdom_sjid == std::vector<int64_t>({0, 99}));
        } else {
            if (specify_domain) {
                REQUIRE(dom_str[0] == dim_infos[1].string_lo);
                REQUIRE(dom_str[1] == dim_infos[1].string_hi);
            } else {
                REQUIRE(dom_str[0] == "");
                REQUIRE(dom_str[1] == "");
            }

            REQUIRE(maxdom_sjid[0] == 0);
            REQUIRE(maxdom_sjid[1] > 2000000000);
        }
        REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));

        sdf->close();

        // Resize
        auto new_shape = int64_t{SOMA_JOINID_RESIZE_DIM_MAX + 1};

        if (!use_current_domain) {
            // Domain is already set. The domain (not current domain but domain)
            // is immutable. All we can do is check for:
            // * throw on write beyond domain
            // * throw on an attempt to resize.
            REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

            sdf = open(OpenMode::write);
            // Array not resizeable if it has not already been sized
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

        } else {
            // Expect throw on write beyond current domain before resize
            REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

            // Check shape after write
            sdf = open(OpenMode::read);
            expect = dim_infos[0].dim_max + 1;
            std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(actual.has_value());
            REQUIRE(actual.value() == expect);
            sdf->close();

            sdf = open(OpenMode::read);
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

            sdf = open(OpenMode::write);
            sdf->resize_soma_joinid_shape(new_shape, "testing");
            sdf->close();

            // Check shape after resize
            sdf = open(OpenMode::read);
            expect = SOMA_JOINID_RESIZE_DIM_MAX + 1;
            actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(actual.has_value());
            REQUIRE(actual.value() == expect);
            sdf->close();

            // Implicitly we expect no throw
            write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX);
        }

        // Check domainish accessors after resize
        sdf->open(OpenMode::read, TimestampRange(0, 2));

        non_empty_domain = sdf->get_non_empty_domain();
        ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            non_empty_domain, "soma_joinid");
        ned_str = ArrowAdapter::get_table_string_column_by_name(
            non_empty_domain, "mystring");

        soma_domain = sdf->get_soma_domain();
        dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            soma_domain, "soma_joinid");
        dom_str = ArrowAdapter::get_table_string_column_by_name(
            soma_domain, "mystring");

        soma_maxdomain = sdf->get_soma_maxdomain();
        maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<
            int64_t>(soma_maxdomain, "soma_joinid");
        maxdom_str = ArrowAdapter::get_table_string_column_by_name(
            soma_maxdomain, "mystring");

        REQUIRE(ned_sjid == std::vector<int64_t>({0, 0}));
        REQUIRE(ned_str == std::vector<std::string>({"", ""}));

        REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));

        if (!use_current_domain) {
            REQUIRE(maxdom_sjid == std::vector<int64_t>({0, 99}));
            REQUIRE(dom_str == std::vector<std::string>({"", ""}));
        } else {
            if (specify_domain) {
                REQUIRE(dom_str[0] == dim_infos[1].string_lo);
                REQUIRE(dom_str[1] == dim_infos[1].string_hi);
            } else {
                REQUIRE(dom_str == std::vector<std::string>({"", ""}));
            }

            REQUIRE(maxdom_sjid[0] == 0);
            REQUIRE(maxdom_sjid[1] > 2000000000);
        }

        REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));

        REQUIRE(ned_str == std::vector<std::string>({"", ""}));

        // Check can_resize_soma_joinid_shape
        std::pair<bool, std::string> check = sdf->can_resize_soma_joinid_shape(
            1, "testing");
        if (!use_current_domain) {
            REQUIRE(check.first == false);
            REQUIRE(
                check.second ==
                "testing: dataframe currently has no domain set: please "
                "upgrade the array.");
        } else {
            // Must fail since this is too small.
            REQUIRE(check.first == false);
            REQUIRE(
                check.second ==
                "testing: new soma_joinid shape 1 < existing shape 99");
            check = sdf->can_resize_soma_joinid_shape(
                SOMA_JOINID_RESIZE_DIM_MAX + 1, "testing");
            REQUIRE(check.first == true);
            REQUIRE(check.second == "");
        }

        sdf->close();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe dim-str-u32 attr-sjid",
    "[SOMADataFrame]") {
    auto use_current_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        auto specify_domain = GENERATE(false, true);
        std::ostringstream section2;
        section2 << "- specify_domain=" << specify_domain;

        std::string suffix1 = use_current_domain ? "true" : "false";
        std::string suffix2 = specify_domain ? "true" : "false";
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-variant-indexed-dataframe-4-" + suffix1 + "-" +
                suffix2);

        std::string string_lo = specify_domain ? "apple" : "";
        std::string string_hi = specify_domain ? "zebra" : "";
        std::vector<helper::DimInfo> dim_infos(
            {str_dim_info(use_current_domain, string_lo, string_hi),
             u32_dim_info(use_current_domain)});
        std::vector<helper::AttrInfo> attr_infos({i64_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::read);

        CurrentDomain current_domain = sdf->get_current_domain_for_test();
        if (!use_current_domain) {
            REQUIRE(current_domain.is_empty());
        } else {
            REQUIRE(!current_domain.is_empty());
            REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
            NDRectangle ndrect = current_domain.ndrectangle();

            std::array<std::string, 2> str_range = ndrect.range<std::string>(
                dim_infos[0].name);
            if (specify_domain) {
                REQUIRE(str_range[0] == dim_infos[0].string_lo);
                REQUIRE(str_range[1] == dim_infos[0].string_hi);
            } else {
                // Can we write empty strings in this range?
                REQUIRE(str_range[0] <= "");
                REQUIRE(str_range[1] >= "");
                // Can we write ASCII values in this range?
                REQUIRE(str_range[0] < " ");
                REQUIRE(str_range[1] > "~");
            }

            std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(
                dim_infos[1].name);
            REQUIRE(u32_range[0] == (uint32_t)0);
            REQUIRE(u32_range[1] == (uint32_t)dim_infos[1].dim_max);
        }

        // Check shape before write
        std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(!actual.has_value());

        // Check domainish accessors before resize
        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<std::string>
            ned_str = ArrowAdapter::get_table_string_column_by_name(
                non_empty_domain, "mystring");

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<std::string>
            dom_str = ArrowAdapter::get_table_string_column_by_name(
                soma_domain, "mystring");

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<std::string>
            maxdom_str = ArrowAdapter::get_table_string_column_by_name(
                soma_maxdomain, "mystring");

        REQUIRE(ned_str == std::vector<std::string>({"", ""}));

        if (!use_current_domain) {
            REQUIRE(dom_str == std::vector<std::string>({"", ""}));
            REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));
        } else {
            if (specify_domain) {
                REQUIRE(dom_str[0] == dim_infos[0].string_lo);
                REQUIRE(dom_str[1] == dim_infos[0].string_hi);
            } else {
                REQUIRE(dom_str == std::vector<std::string>({"", ""}));
            }
            REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));
        }

        sdf->close();

        REQUIRE(sdf->nnz() == 0);

        // Write
        write_sjid_u32_str_data_from(0);

        REQUIRE(sdf->nnz() == 2);
        write_sjid_u32_str_data_from(8);
        // soma_joinid is not a dim here and so the second write is an overwrite
        // of the first here
        REQUIRE(sdf->nnz() == 2);

        // Check shape after write
        sdf = open(OpenMode::read);
        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(!actual.has_value());
        sdf->close();

        // Resize
        auto new_shape = int64_t{SOMA_JOINID_RESIZE_DIM_MAX + 1};

        if (!use_current_domain) {
            // Domain is already set. The domain (not current domain but domain)
            // is immutable. All we can do is check for:
            // * throw on write beyond domain -- except here, soma_joinid is not
            //   a dim, so no throw
            // * throw on an attempt to resize.

            sdf = open(OpenMode::write);
            // Array not resizeable if it has not already been sized
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

        } else {
            // Check shape after write
            sdf = open(OpenMode::read);
            std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(!actual.has_value());
            sdf->close();

            sdf = open(OpenMode::read);
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

            sdf = open(OpenMode::write);
            sdf->resize_soma_joinid_shape(new_shape, "testing");
            sdf->close();

            // Check shape after resize -- noting soma_joinid is not a dim here
            sdf = open(OpenMode::read);
            actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(!actual.has_value());
            sdf->close();

            // Implicitly we expect no throw
            write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX);
        }

        // Check domainish accessors after resize
        sdf->open(OpenMode::read, TimestampRange(0, 2));

        non_empty_domain = sdf->get_non_empty_domain();
        ned_str = ArrowAdapter::get_table_string_column_by_name(
            non_empty_domain, "mystring");

        soma_domain = sdf->get_soma_domain();
        dom_str = ArrowAdapter::get_table_string_column_by_name(
            soma_domain, "mystring");

        soma_maxdomain = sdf->get_soma_maxdomain();
        maxdom_str = ArrowAdapter::get_table_string_column_by_name(
            soma_maxdomain, "mystring");

        REQUIRE(ned_str == std::vector<std::string>({"", ""}));

        if (!use_current_domain) {
            REQUIRE(dom_str == std::vector<std::string>({"", ""}));
            REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));
        } else {
            if (specify_domain) {
                REQUIRE(dom_str[0] == dim_infos[0].string_lo);
                REQUIRE(dom_str[1] == dim_infos[0].string_hi);
            } else {
                REQUIRE(dom_str == std::vector<std::string>({"", ""}));
            }
            REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));
        }

        // Check can_resize_soma_joinid_shape
        std::pair<bool, std::string> check = sdf->can_resize_soma_joinid_shape(
            0, "testing");
        if (!use_current_domain) {
            REQUIRE(check.first == false);
            REQUIRE(
                check.second ==
                "testing: dataframe currently has no domain set: please "
                "upgrade the array.");
        } else {
            // Must pass since soma_joinid isn't a dim in this case.
            REQUIRE(check.first == true);
            REQUIRE(check.second == "");
        }

        sdf->close();
    }
}
