/**
 * @file   unit_sdf.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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

    std::string i64_arrow_format = ArrowAdapter::tdb_to_arrow_type(i64_datatype);
    std::string u32_arrow_format = ArrowAdapter::tdb_to_arrow_type(u32_datatype);
    std::string attr_1_arrow_format = ArrowAdapter::tdb_to_arrow_type(str_datatype);

    helper::DimInfo i64_dim_info() {
        return helper::DimInfo(
            {.name = i64_name,
             .tiledb_datatype = i64_datatype,
             .dim_max = i64_dim_max,
             .string_lo = "N/A",
             .string_hi = "N/A"});
    }
    helper::DimInfo u32_dim_info() {
        return helper::DimInfo(
            {.name = u32_name,
             .tiledb_datatype = u32_datatype,
             .dim_max = u32_dim_max,
             .string_lo = "N/A",
             .string_hi = "N/A"});
    }
    helper::DimInfo str_dim_info(std::string string_lo, std::string string_hi) {
        return helper::DimInfo(
            {.name = str_name,
             .tiledb_datatype = str_datatype,
             .dim_max = str_dim_max,
             .string_lo = string_lo,
             .string_hi = string_hi});
    }

    helper::AttrInfo i64_attr_info(std::string name = i64_name) {
        return helper::AttrInfo({.name = name, .tiledb_datatype = i64_datatype});
    }
    helper::AttrInfo u32_attr_info() {
        return helper::AttrInfo({.name = u32_name, .tiledb_datatype = u32_datatype});
    }
    helper::AttrInfo str_attr_info() {
        return helper::AttrInfo({.name = str_name, .tiledb_datatype = str_datatype});
    }

    //  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Helper methods for create/open/write/etc.

    void create(const std::vector<helper::DimInfo>& dim_infos, const std::vector<helper::AttrInfo>& attr_infos) {
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(dim_infos, attr_infos);
        SOMADataFrame::create(uri_, schema, index_columns, ctx_);
    }

    void create(
        const std::vector<helper::DimInfo>& dim_infos,
        const std::vector<helper::AttrInfo>& attr_infos,
        const PlatformConfig& platform_config,
        std::optional<TimestampRange> timestamp_range = std::nullopt) {
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(dim_infos, attr_infos);
        SOMADataFrame::create(uri_, schema, index_columns, ctx_, platform_config, timestamp_range);
    }

    std::unique_ptr<SOMADataFrame> open(OpenMode mode, std::optional<TimestampRange> timestamp_range = std::nullopt) {
        return SOMADataFrame::open(uri_, mode, ctx_, timestamp_range);
    }

    void write_sjid_u32_str_data_from(int64_t sjid_base) {
        auto sdf = SOMADataFrame::open(uri_, OpenMode::soma_write, ctx_);

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

        auto mq = sdf->create_managed_query();
        mq.setup_write_column(i64_name, i64_data.size(), i64_data.data(), (uint64_t*)nullptr);
        mq.setup_write_column(str_name, strings.size(), char_data.data(), char_offsets.data());
        mq.setup_write_column(u32_name, u32_data.size(), u32_data.data(), (uint64_t*)nullptr);
        mq.submit_write();

        sdf->close();
    }
};

TEST_CASE_METHOD(VariouslyIndexedDataFrameFixture, "SOMADataFrame: basic", "[SOMADataFrame]") {
    set_up(std::make_shared<SOMAContext>(), "mem://unit-test-dataframe-basic");

    std::vector<helper::DimInfo> dim_infos({i64_dim_info()});
    std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

    REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

    create(dim_infos, attr_infos);

    REQUIRE(SOMADataFrame::exists(uri_, ctx_));
    REQUIRE(!SOMASparseNDArray::exists(uri_, ctx_));
    REQUIRE(!SOMADenseNDArray::exists(uri_, ctx_));

    auto sdf = open(OpenMode::soma_read);
    REQUIRE(sdf->uri() == uri_);
    REQUIRE(sdf->ctx() == ctx_);
    REQUIRE(sdf->type() == "SOMADataFrame");
    std::vector<std::string> expected_index_column_names = {dim_infos[0].name};
    REQUIRE(sdf->index_column_names() == expected_index_column_names);
    REQUIRE(sdf->nnz() == 0);
    sdf->close();

    std::vector<int64_t> d0(10);
    for (int j = 0; j < 10; j++)
        d0[j] = j;
    std::vector<uint32_t> a0(10, 1);

    // A write in read mode should fail
    {
        sdf = open(OpenMode::soma_read);
        auto mq = sdf->create_managed_query();
        REQUIRE_THROWS(mq.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr));
        REQUIRE_THROWS(mq.setup_write_column(attr_infos[0].name, a0.size(), a0.data(), (uint64_t*)nullptr));
        REQUIRE_THROWS(mq.submit_write());
        sdf->close();
    }

    {
        sdf = open(OpenMode::soma_write);
        auto mq = sdf->create_managed_query();
        mq.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column(attr_infos[0].name, a0.size(), a0.data(), (uint64_t*)nullptr);
        mq.submit_write();
        sdf->close();
    }

    {
        sdf = open(OpenMode::soma_read);
        auto mq = sdf->create_managed_query();
        while (auto batch = mq.read_next()) {
            auto arrbuf = batch.value();
            auto d0span = arrbuf->at(dim_infos[0].name)->data<int64_t>();
            auto a0span = arrbuf->at(attr_infos[0].name)->data<uint32_t>();
            REQUIRE(d0 == std::vector<int64_t>(d0span.begin(), d0span.end()));
            REQUIRE(a0 == std::vector<uint32_t>(a0span.begin(), a0span.end()));
        }
        sdf->close();
    }

    auto soma_object = SOMAObject::open(uri_, OpenMode::soma_read, ctx_);
    REQUIRE(soma_object->uri() == uri_);
    REQUIRE(soma_object->type() == "SOMADataFrame");
    soma_object->close();
}

TEST_CASE_METHOD(VariouslyIndexedDataFrameFixture, "SOMADataFrame: platform_config", "[SOMADataFrame]") {
    std::pair<std::string, tiledb_filter_type_t> filter = GENERATE(
        std::make_pair(R"({"name": "GZIP", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_GZIP),
        std::make_pair(R"({"name": "ZSTD", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_ZSTD),
        std::make_pair(R"({"name": "LZ4", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_LZ4),
        std::make_pair(R"({"name": "BZIP2", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_BZIP2),
        std::make_pair(R"({"name": "RLE", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_RLE),
        std::make_pair(R"({"name": "DICTIONARY_ENCODING", "COMPRESSION_LEVEL": 3})", TILEDB_FILTER_DICTIONARY),
        std::make_pair(
            R"({"name": "BIT_WIDTH_REDUCTION", "BIT_WIDTH_MAX_WINDOW": 3})", TILEDB_FILTER_BIT_WIDTH_REDUCTION),
        std::make_pair(R"({"name": "POSITIVE_DELTA", "POSITIVE_DELTA_MAX_WINDOW": 3})", TILEDB_FILTER_POSITIVE_DELTA),
        std::make_pair(
            R"({"name": "DELTA", "COMPRESSION_LEVEL": 3, "COMPRESSION_REINTERPRET_DATATYPE": 1})", TILEDB_FILTER_DELTA),
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

    SECTION("- filter=" + filter.first) {
        set_up(std::make_shared<SOMAContext>(), "mem://unit-test-dataframe-platform-config");

        PlatformConfig platform_config;
        platform_config.cell_order = "hilbert";
        platform_config.dataframe_dim_zstd_level = 6;
        platform_config.offsets_filters = R"([)" + filter.first + R"(])";
        platform_config.validity_filters = R"([)" + filter.first + R"(])";
        if (filter.second != TILEDB_FILTER_WEBP) {
            platform_config.attrs = R"({"a0": {"filters":[)" + filter.first + R"(]}})";
        }

        std::vector<helper::DimInfo> dim_infos({i64_dim_info()});
        std::vector<helper::AttrInfo> attr_infos({i64_attr_info("a0")});

        REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

        create(dim_infos, attr_infos, platform_config);

        auto sdf = open(OpenMode::soma_read);
        auto sch = sdf->tiledb_schema();
        REQUIRE(sch->offsets_filter_list().filter(0).filter_type() == filter.second);

        REQUIRE(sch->validity_filter_list().filter(0).filter_type() == filter.second);

        auto dim_filter = sch->domain().dimension(dim_infos[0].name).filter_list().filter(0);
        REQUIRE(dim_filter.filter_type() == TILEDB_FILTER_ZSTD);
        REQUIRE(dim_filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL) == 6);

        if (filter.second != TILEDB_FILTER_WEBP) {
            REQUIRE(sch->attribute(attr_infos[0].name).filter_list().filter(0).filter_type() == filter.second);
        }

        auto config_options = sdf->schema_config_options();
        REQUIRE(config_options.capacity == 100000);
        REQUIRE(config_options.allows_duplicates == false);
        REQUIRE(config_options.tile_order == "row-major");
        REQUIRE(config_options.cell_order == "hilbert");

        REQUIRE(json::parse(config_options.offsets_filters)[0].at("name") == Filter::to_str(filter.second));
        REQUIRE(json::parse(config_options.validity_filters)[0].at("name") == Filter::to_str(filter.second));
        if (filter.second != TILEDB_FILTER_WEBP) {
            REQUIRE(json::parse(config_options.attrs)["a0"]["filters"][0].at("name") == Filter::to_str(filter.second));
        }
        REQUIRE(
            json::parse(config_options.dims)["soma_joinid"]["filters"][0].at("name") ==
            Filter::to_str(TILEDB_FILTER_ZSTD));

        sdf->close();
    }
}

TEST_CASE_METHOD(VariouslyIndexedDataFrameFixture, "SOMADataFrame: metadata", "[SOMADataFrame]") {
    set_up(std::make_shared<SOMAContext>(), "mem://unit-test-collection");

    std::vector<helper::DimInfo> dim_infos({i64_dim_info()});
    std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

    REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

    create(dim_infos, attr_infos, PlatformConfig(), TimestampRange(0, 1));

    auto sdf = open(OpenMode::soma_write, TimestampRange(0, 2));

    int32_t val = 100;
    sdf->set_metadata("md", TILEDB_INT32, 1, &val);
    sdf->close();

    // Read metadata
    sdf->open(OpenMode::soma_read, TimestampRange(0, 2));
    REQUIRE(sdf->metadata_num() == 3);
    REQUIRE(sdf->has_metadata("soma_object_type"));
    REQUIRE(sdf->has_metadata("soma_encoding_version"));
    REQUIRE(sdf->has_metadata("md"));
    auto mdval = sdf->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    sdf->close();

    // md should not be available at (0, 1)
    sdf->open(OpenMode::soma_read, TimestampRange(0, 1));
    REQUIRE(sdf->metadata_num() == 2);
    REQUIRE(sdf->has_metadata("soma_object_type"));
    REQUIRE(sdf->has_metadata("soma_encoding_version"));
    REQUIRE(!sdf->has_metadata("md"));
    sdf->close();

    // Metadata should also be retrievable in write mode
    sdf->open(OpenMode::soma_write);
    REQUIRE(sdf->metadata_num() == 3);
    REQUIRE(sdf->has_metadata("soma_object_type"));
    REQUIRE(sdf->has_metadata("soma_encoding_version"));
    REQUIRE(sdf->has_metadata("md"));
    mdval = sdf->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    sdf->delete_metadata("md");
    mdval = sdf->get_metadata("md");
    REQUIRE(!mdval.has_value());
    sdf->close();

    // Confirm delete in read mode
    sdf->open(OpenMode::soma_read);
    REQUIRE(!sdf->has_metadata("md"));
    REQUIRE(sdf->metadata_num() == 2);
}

TEST_CASE_METHOD(VariouslyIndexedDataFrameFixture, "SOMADataFrame: bounds-checking", "[SOMADataFrame]") {
    int old_max = SOMA_JOINID_DIM_MAX;
    int new_max = SOMA_JOINID_RESIZE_DIM_MAX;

    set_up(std::make_shared<SOMAContext>(), "mem://unit-test-bounds-checking");

    std::vector<helper::DimInfo> dim_infos({i64_dim_info()});
    std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

    REQUIRE(!SOMADataFrame::exists(uri_, ctx_));

    create(dim_infos, attr_infos);
    std::vector<int64_t> d0({old_max + 1, old_max + 2});
    std::vector<double> a0({1.5, 2.5});

    {
        auto sdf = open(OpenMode::soma_write);
        auto mq = sdf->create_managed_query();

        mq.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column(attr_infos[0].name, a0.size(), a0.data(), (uint64_t*)nullptr);
        // Writing outside the current domain should fail
        REQUIRE_THROWS(mq.submit_write());
        sdf->close();

        sdf = open(OpenMode::soma_write);
        sdf->resize_soma_joinid_shape(int64_t{new_max}, "testing");
        sdf->close();
    }

    {
        auto sdf = open(OpenMode::soma_write);
        auto mq = sdf->create_managed_query();
        mq.setup_write_column(dim_infos[0].name, d0.size(), d0.data(), (uint64_t*)nullptr);
        mq.setup_write_column(attr_infos[0].name, a0.size(), a0.data(), (uint64_t*)nullptr);
        // Writing after resize should succeed
        mq.submit_write();

        sdf->close();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: standard-indexed dataframe dim-sjid attr-str-u32",
    "[SOMADataFrame]") {
    // We have these:
    // * upgrade_domain requires the user to specify values for all index
    //   columns. This is in the spec.
    // * resize_soma_joinid_shape allows the user to specify only the
    //   desired soma_joinid shape. This is crucial for experiment-level
    //   resize as an internal method at the Python level.
    // Both need testing. Each one adds a shape where there wasn't one
    // before. So we need to test one or the other on a given run.
    auto test_upgrade_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- test_upgrade_domain=" << test_upgrade_domain;
    SECTION(section.str()) {
        std::string suffix = test_upgrade_domain ? "true" : "false";

        set_up(std::make_shared<SOMAContext>(), "mem://unit-test-variant-indexed-dataframe-1-" + suffix);

        std::vector<helper::DimInfo> dim_infos({i64_dim_info()});
        std::vector<helper::AttrInfo> attr_infos({str_attr_info(), u32_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::soma_read);

        CurrentDomain current_domain = sdf->get_current_domain_for_test();
        REQUIRE(!current_domain.is_empty());
        REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
        NDRectangle ndrect = current_domain.ndrectangle();

        std::array<int64_t, 2> i64_range = ndrect.range<int64_t>(dim_infos[0].name);
        REQUIRE(i64_range[0] == (int64_t)0);
        REQUIRE(i64_range[1] == (int64_t)dim_infos[0].dim_max);

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
        sdf->open(OpenMode::soma_read);

        REQUIRE(sdf->nnz() == 2);

        expect = dim_infos[0].dim_max + 1;
        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);

        // Check domainish accessors before resize
        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<int64_t> ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            non_empty_domain, "soma_joinid");

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<int64_t> dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            soma_domain, "soma_joinid");

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<int64_t> maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            soma_maxdomain, "soma_joinid");

        REQUIRE(ned_sjid == std::vector<int64_t>({1, 2}));

        REQUIRE(dom_sjid == std::vector<int64_t>({0, SOMA_JOINID_DIM_MAX}));

        REQUIRE(maxdom_sjid.size() == 2);
        REQUIRE(maxdom_sjid[0] == 0);
        REQUIRE(maxdom_sjid[1] > 2000000000);
        sdf->close();

        REQUIRE(sdf->nnz() == 2);
        write_sjid_u32_str_data_from(8);
        REQUIRE(sdf->nnz() == 4);

        sdf->open(OpenMode::soma_read);

        // Check can_upgrade_domain
        // std::unique_ptr<ArrowSchema>
        //     domain_schema = create_index_cols_info_schema(dim_infos);
        // auto domain_array = ArrowAdapter::make_arrow_array_parent(
        //     dim_infos.size());
        // // OK since there currently is no shape set:
        // domain_array->children[0] = ArrowAdapter::make_arrow_array_child(
        //     std::vector<int64_t>({0, 0}));
        // auto domain_table = ArrowTable(
        //     std::move(domain_array), std::move(domain_schema));
        StatusAndReason check = sdf->can_upgrade_soma_joinid_shape(1, "testing");
        // Must fail since this is too small.
        REQUIRE(check.first == false);
        REQUIRE(check.second == "testing: dataframe already has its domain set.");

        // Check can_upgrade_soma_joinid_shape
        check = sdf->can_upgrade_soma_joinid_shape(1, "testing");
        // Must fail since this is too small.
        REQUIRE(check.first == false);
        REQUIRE(check.second == "testing: dataframe already has its domain set.");

        sdf->close();

        // Resize
        auto new_shape = int64_t{SOMA_JOINID_RESIZE_DIM_MAX + 1};

        // Expect throw on write beyond current domain before resize
        REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

        // Check shape after write
        sdf = open(OpenMode::soma_read);
        expect = dim_infos[0].dim_max + 1;

        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);
        sdf->close();

        // Apply the domain change
        if (test_upgrade_domain) {
            managed_unique_ptr<ArrowSchema> domain_schema = create_index_cols_info_schema(dim_infos);
            auto domain_array = ArrowAdapter::make_arrow_array_parent(dim_infos.size());
            domain_array->children[0] = ArrowAdapter::make_arrow_array_child(std::vector<int64_t>({0, new_shape - 1}));
            auto domain_table = ArrowTable(std::move(domain_array), std::move(domain_schema));

            // Not open for write
            sdf = open(OpenMode::soma_read);
            REQUIRE_THROWS(sdf->change_domain(domain_table, "testing"));
            sdf->close();

            // Open for write
            sdf = open(OpenMode::soma_write);
            sdf->change_domain(domain_table, "testing");
            sdf->close();
        } else {
            // Not open for write
            sdf = open(OpenMode::soma_read);
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

            // Open for write
            sdf = open(OpenMode::soma_write);
            sdf->resize_soma_joinid_shape(new_shape, "testing");
            sdf->close();
        }

        // Check shape after resize
        sdf = open(OpenMode::soma_read);
        expect = SOMA_JOINID_RESIZE_DIM_MAX + 1;
        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);
        sdf->close();

        // Implicitly we expect no throw
        write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX);

        // Check domainish accessors after resize
        sdf->open(OpenMode::soma_read);

        non_empty_domain = sdf->get_non_empty_domain();
        ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(non_empty_domain, "soma_joinid");

        soma_domain = sdf->get_soma_domain();
        dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(soma_domain, "soma_joinid");

        soma_maxdomain = sdf->get_soma_maxdomain();
        maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(soma_maxdomain, "soma_joinid");

        REQUIRE(ned_sjid == std::vector<int64_t>({1, 101}));
        REQUIRE(dom_sjid == std::vector<int64_t>({0, SOMA_JOINID_RESIZE_DIM_MAX}));
        REQUIRE(maxdom_sjid.size() == 2);
        REQUIRE(maxdom_sjid[0] == 0);
        REQUIRE(maxdom_sjid[1] > 2000000000);

        // Check can_resize_soma_joinid_shape
        check = sdf->can_resize_soma_joinid_shape(1, "testing");
        // Must fail since this is too small.
        REQUIRE(check.first == false);
        REQUIRE(
            check.second ==
            "testing: new soma_joinid shape 1 < existing shape "
            "200");
        check = sdf->can_resize_soma_joinid_shape(SOMA_JOINID_RESIZE_DIM_MAX + 1, "testing");
        REQUIRE(check.first == true);
        REQUIRE(check.second == "");

        sdf->close();

        // Check can_upgrade_domain
        sdf->open(OpenMode::soma_read);
        // The dataframe now has a shape
        check = sdf->can_upgrade_soma_joinid_shape(1, "testing");
        // Must fail since this is too small.
        REQUIRE(check.first == false);
        REQUIRE(check.second == "testing: dataframe already has its domain set.");
        sdf->close();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe dim-u32-sjid attr-str",
    "[SOMADataFrame]") {
    // We have these:
    // * upgrade_domain requires the user to specify values for all
    // index
    //   columns. This is in the spec.
    // * resize_soma_joinid_shape allows the user to specify only the
    //   desired soma_joinid shape. This is crucial for experiment-level
    //   resize as an internal method at the Python level.
    // Both need testing. Each one adds a shape where there wasn't one
    // before. So we need to test one or the other on a given run.
    auto test_upgrade_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- test_upgrade_domain=" << test_upgrade_domain;
    SECTION(section.str()) {
        std::string suffix = test_upgrade_domain ? "true" : "false";
        set_up(std::make_shared<SOMAContext>(), "mem://unit-test-variant-indexed-dataframe-2-" + suffix);

        std::vector<helper::DimInfo> dim_infos({u32_dim_info(), i64_dim_info()});
        std::vector<helper::AttrInfo> attr_infos({str_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::soma_read);

        CurrentDomain current_domain = sdf->get_current_domain_for_test();
        REQUIRE(!current_domain.is_empty());
        REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
        NDRectangle ndrect = current_domain.ndrectangle();

        std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(dim_infos[0].name);
        REQUIRE(u32_range[0] == (uint32_t)0);
        REQUIRE(u32_range[1] == (uint32_t)dim_infos[0].dim_max);

        std::array<int64_t, 2> i64_range = ndrect.range<int64_t>(dim_infos[1].name);
        REQUIRE(i64_range[0] == (int64_t)0);
        REQUIRE(i64_range[1] == (int64_t)dim_infos[1].dim_max);

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
        sdf->open(OpenMode::soma_read);

        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<int64_t> ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            non_empty_domain, "soma_joinid");
        std::vector<uint32_t> ned_u32 = ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(
            non_empty_domain, "myuint32");

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<int64_t> dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            soma_domain, "soma_joinid");
        std::vector<uint32_t> dom_u32 = ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(
            soma_domain, "myuint32");

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<int64_t> maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
            soma_maxdomain, "soma_joinid");
        std::vector<uint32_t> maxdom_u32 = ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(
            soma_maxdomain, "myuint32");

        REQUIRE(ned_sjid == std::vector<int64_t>({1, 10}));
        REQUIRE(ned_u32 == std::vector<uint32_t>({1234, 5678}));

        REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));
        REQUIRE(dom_u32 == std::vector<uint32_t>({0, 9999}));

        REQUIRE(maxdom_sjid.size() == 2);
        REQUIRE(maxdom_u32.size() == 2);

        REQUIRE(maxdom_u32[0] == 0);
        REQUIRE(maxdom_u32[1] > 2000000000);

        sdf->close();

        // Check shape after write
        sdf = open(OpenMode::soma_read);
        expect = dim_infos[1].dim_max + 1;
        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);

        // Check can_upgrade_soma_joinid_shape
        StatusAndReason check = sdf->can_upgrade_soma_joinid_shape(1, "testing");
        // Must fail since this is too small.
        REQUIRE(check.first == false);
        REQUIRE(check.second == "testing: dataframe already has its domain set.");

        // Check can_upgrade_domain
        check = sdf->can_upgrade_soma_joinid_shape(1, "testing");
        // Must fail since this is too small.
        REQUIRE(check.first == false);
        REQUIRE(check.second == "testing: dataframe already has its domain set.");

        sdf->close();

        // Resize
        auto new_shape = int64_t{SOMA_JOINID_RESIZE_DIM_MAX + 1};
        uint32_t new_u32_dim_max = (uint32_t)u32_dim_max * 2 + 1;

        // Expect throw on write beyond current domain before resize
        REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

        // Check shape after write
        sdf = open(OpenMode::soma_read);
        expect = dim_infos[1].dim_max + 1;
        actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(actual.has_value());
        REQUIRE(actual.value() == expect);
        sdf->close();

        // Apply the domain change
        if (test_upgrade_domain) {
            managed_unique_ptr<ArrowSchema> domain_schema = create_index_cols_info_schema(dim_infos);
            auto domain_array = ArrowAdapter::make_arrow_array_parent(dim_infos.size());
            domain_array->children[0] = ArrowAdapter::make_arrow_array_child(
                std::vector<uint32_t>({0, new_u32_dim_max}));
            domain_array->children[1] = ArrowAdapter::make_arrow_array_child(std::vector<int64_t>({0, new_shape - 1}));
            auto domain_table = ArrowTable(std::move(domain_array), std::move(domain_schema));

            // Not open for write
            sdf = open(OpenMode::soma_read);
            REQUIRE_THROWS(sdf->change_domain(domain_table, "testing"));
            sdf->close();

            // Open for write
            sdf = open(OpenMode::soma_write);
            sdf->change_domain(domain_table, "testing");
            sdf->close();
        } else {
            // Not open for write
            sdf = open(OpenMode::soma_read);
            REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
            sdf->close();

            // Open for write
            sdf = open(OpenMode::soma_write);
            sdf->resize_soma_joinid_shape(new_shape, "testing");
            sdf->close();
        }
        sdf->close();

        // Implicitly we expect no throw
        write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX);

        // Check domainish accessors after resize
        sdf->open(OpenMode::soma_read);

        non_empty_domain = sdf->get_non_empty_domain();
        ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(non_empty_domain, "soma_joinid");
        ned_u32 = ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(non_empty_domain, "myuint32");

        soma_domain = sdf->get_soma_domain();
        dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(soma_domain, "soma_joinid");
        dom_u32 = ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(soma_domain, "myuint32");

        soma_maxdomain = sdf->get_soma_maxdomain();
        maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(soma_maxdomain, "soma_joinid");
        maxdom_u32 = ArrowAdapter::get_table_non_string_column_by_name<uint32_t>(soma_maxdomain, "myuint32");

        REQUIRE(ned_sjid == std::vector<int64_t>({1, 101}));
        REQUIRE(ned_u32 == std::vector<uint32_t>({1234, 5678}));

        REQUIRE(dom_sjid == std::vector<int64_t>({0, 199}));
        if (test_upgrade_domain) {
            REQUIRE(dom_u32 == std::vector<uint32_t>({0, 19999}));
        } else {
            REQUIRE(dom_u32 == std::vector<uint32_t>({0, 9999}));
        }

        REQUIRE(maxdom_sjid.size() == 2);
        REQUIRE(maxdom_sjid[0] == 0);
        REQUIRE(maxdom_sjid[1] > 2000000000);

        REQUIRE(maxdom_u32.size() == 2);
        REQUIRE(maxdom_u32[0] == 0);
        REQUIRE(maxdom_u32[1] > 2000000000);

        // Check can_resize_soma_joinid_shape
        check = sdf->can_resize_soma_joinid_shape(1, "testing");
        // Must fail since this is too small.
        REQUIRE(check.first == false);
        REQUIRE(
            check.second ==
            "testing: new soma_joinid shape 1 < existing shape "
            "200");
        check = sdf->can_resize_soma_joinid_shape(SOMA_JOINID_RESIZE_DIM_MAX + 1, "testing");
        REQUIRE(check.first == true);
        REQUIRE(check.second == "");

        sdf->close();
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe dim-sjid-str attr-u32",
    "[SOMADataFrame]") {
    auto specify_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- specify_domain=" << specify_domain;
    SECTION(section.str()) {
        auto test_upgrade_domain = GENERATE(false, true);
        std::ostringstream section3;
        section << "- test_upgrade_domain=" << test_upgrade_domain;
        SECTION(section3.str()) {
            std::string suffix1 = specify_domain ? "true" : "false";
            std::string suffix2 = test_upgrade_domain ? "true" : "false";
            set_up(
                std::make_shared<SOMAContext>(),
                "mem://unit-test-variant-indexed-dataframe-3-" + suffix1 + "-" + suffix2);

            std::string string_lo = "";
            std::string string_hi = "";
            std::vector<helper::DimInfo> dim_infos({i64_dim_info(), str_dim_info(string_lo, string_hi)});
            std::vector<helper::AttrInfo> attr_infos({u32_attr_info()});

            // Create
            create(dim_infos, attr_infos);

            // Check current domain
            auto sdf = open(OpenMode::soma_read);

            CurrentDomain current_domain = sdf->get_current_domain_for_test();
            REQUIRE(!current_domain.is_empty());
            REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
            NDRectangle ndrect = current_domain.ndrectangle();

            std::array<int64_t, 2> i64_range = ndrect.range<int64_t>(dim_infos[0].name);
            REQUIRE(i64_range[0] == (int64_t)0);
            REQUIRE(i64_range[1] == (int64_t)dim_infos[0].dim_max);

            std::array<std::string, 2> str_range = ndrect.range<std::string>(dim_infos[1].name);

            // Can we write empty strings in this range?
            REQUIRE(str_range[0] <= "");
            REQUIRE(str_range[1] >= "");
            // Can we write ASCII values in this range?
            REQUIRE(str_range[0] < " ");
            REQUIRE(str_range[1] > "~");

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
            sdf = open(OpenMode::soma_read);
            expect = dim_infos[0].dim_max + 1;
            actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(actual.has_value());
            REQUIRE(actual.value() == expect);

            // Check domainish accessors before resize
            ArrowTable non_empty_domain = sdf->get_non_empty_domain();
            std::vector<int64_t> ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                non_empty_domain, "soma_joinid");
            std::vector<std::string> ned_str = ArrowAdapter::get_table_string_column_by_name(
                non_empty_domain, "mystring");

            ArrowTable soma_domain = sdf->get_soma_domain();
            std::vector<int64_t> dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                soma_domain, "soma_joinid");
            std::vector<std::string> dom_str = ArrowAdapter::get_table_string_column_by_name(soma_domain, "mystring");

            ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
            std::vector<int64_t> maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
                soma_maxdomain, "soma_joinid");
            std::vector<std::string> maxdom_str = ArrowAdapter::get_table_string_column_by_name(
                soma_maxdomain, "mystring");

            REQUIRE(ned_sjid == std::vector<int64_t>({1, 10}));
            REQUIRE(ned_str == std::vector<std::string>({"apple", "bat"}));

            REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));

            if (specify_domain) {
                REQUIRE(dom_str[0] == dim_infos[1].string_lo);
                REQUIRE(dom_str[1] == dim_infos[1].string_hi);
            } else {
                REQUIRE(dom_str[0] == "");
                REQUIRE(dom_str[1] == "");
            }

            REQUIRE(maxdom_sjid[0] == 0);
            REQUIRE(maxdom_sjid[1] > 2000000000);
            REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));

            sdf->close();

            // Check can_upgrade_domain
            sdf = open(OpenMode::soma_read);
            StatusAndReason check = sdf->can_upgrade_soma_joinid_shape(1, "testing");
            // Must fail since this is too small.
            REQUIRE(check.first == false);
            REQUIRE(check.second == "testing: dataframe already has its domain set.");

            sdf->close();

            // Resize

            auto new_shape = int64_t{SOMA_JOINID_RESIZE_DIM_MAX + 1};

            // Expect throw on write beyond current domain before
            // resize
            REQUIRE_THROWS(write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX));

            // Check shape after write
            sdf = open(OpenMode::soma_read);
            expect = dim_infos[0].dim_max + 1;
            actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(actual.has_value());
            REQUIRE(actual.value() == expect);
            sdf->close();

            // Apply the domain change
            if (test_upgrade_domain) {
                managed_unique_ptr<ArrowSchema> domain_schema = create_index_cols_info_schema(dim_infos);
                auto domain_array = ArrowAdapter::make_arrow_array_parent(dim_infos.size());
                domain_array->children[0] = ArrowAdapter::make_arrow_array_child(
                    std::vector<int64_t>({0, new_shape - 1}));
                domain_array->children[1] = ArrowAdapter::make_arrow_array_child_string(
                    std::vector<std::string>({"", ""}));
                auto domain_table = ArrowTable(std::move(domain_array), std::move(domain_schema));

                // Not open for write
                sdf = open(OpenMode::soma_read);
                REQUIRE_THROWS(sdf->change_domain(domain_table, "testing"));
                sdf->close();

                // Open for write
                sdf = open(OpenMode::soma_write);
                sdf->change_domain(domain_table, "testing");
                sdf->close();
            } else {
                // Not open for write
                sdf = open(OpenMode::soma_read);
                REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
                sdf->close();

                // Open for write
                sdf = open(OpenMode::soma_write);
                sdf->resize_soma_joinid_shape(new_shape, "testing");

                sdf->close();
            }

            sdf->open(OpenMode::soma_write);
            // Implicitly we expect no throw
            write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX);
            sdf->close();

            // Check domainish accessors after resize
            sdf->open(OpenMode::soma_read, TimestampRange(0, 2));

            non_empty_domain = sdf->get_non_empty_domain();
            ned_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(non_empty_domain, "soma_joinid");
            ned_str = ArrowAdapter::get_table_string_column_by_name(non_empty_domain, "mystring");

            soma_domain = sdf->get_soma_domain();
            dom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(soma_domain, "soma_joinid");
            dom_str = ArrowAdapter::get_table_string_column_by_name(soma_domain, "mystring");

            soma_maxdomain = sdf->get_soma_maxdomain();
            maxdom_sjid = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(soma_maxdomain, "soma_joinid");
            maxdom_str = ArrowAdapter::get_table_string_column_by_name(soma_maxdomain, "mystring");

            REQUIRE(ned_sjid == std::vector<int64_t>({0, 0}));
            REQUIRE(ned_str == std::vector<std::string>({"", ""}));

            REQUIRE(dom_sjid == std::vector<int64_t>({0, 99}));

            if (specify_domain) {
                REQUIRE(dom_str[0] == dim_infos[1].string_lo);
                REQUIRE(dom_str[1] == dim_infos[1].string_hi);
            } else {
                REQUIRE(dom_str == std::vector<std::string>({"", ""}));
            }

            REQUIRE(maxdom_sjid[0] == 0);
            REQUIRE(maxdom_sjid[1] > 2000000000);

            REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));

            REQUIRE(ned_str == std::vector<std::string>({"", ""}));

            // Check can_resize_soma_joinid_shape
            check = sdf->can_resize_soma_joinid_shape(1, "testing");
            // Must fail since this is too small.
            REQUIRE(check.first == false);
            REQUIRE(
                check.second ==
                "testing: new soma_joinid shape 1 < existing shape "
                "100");
            check = sdf->can_resize_soma_joinid_shape(SOMA_JOINID_RESIZE_DIM_MAX + 1, "testing");
            REQUIRE(check.first == true);
            REQUIRE(check.second == "");

            sdf->close();
        }
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMADataFrame: variant-indexed dataframe dim-str-u32 attr-sjid",
    "[SOMADataFrame]") {
    auto specify_domain = GENERATE(false, true);
    std::ostringstream section;
    section << "- specify_domain=" << specify_domain;
    SECTION(section.str()) {
        auto test_upgrade_domain = GENERATE(false, true);
        std::ostringstream section3;
        section << "- test_upgrade_domain=" << test_upgrade_domain;
        SECTION(section3.str()) {
            std::string suffix1 = specify_domain ? "true" : "false";
            std::string suffix2 = test_upgrade_domain ? "true" : "false";
            set_up(
                std::make_shared<SOMAContext>(),
                "mem://unit-test-variant-indexed-dataframe-4-" + suffix1 + "-" + suffix2);

            std::string string_lo = "";
            std::string string_hi = "";
            std::vector<helper::DimInfo> dim_infos({str_dim_info(string_lo, string_hi), u32_dim_info()});
            std::vector<helper::AttrInfo> attr_infos({i64_attr_info()});

            // Create
            create(dim_infos, attr_infos);

            // Check current domain
            auto sdf = open(OpenMode::soma_read);

            CurrentDomain current_domain = sdf->get_current_domain_for_test();
            REQUIRE(!current_domain.is_empty());
            REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
            NDRectangle ndrect = current_domain.ndrectangle();

            std::array<std::string, 2> str_range = ndrect.range<std::string>(dim_infos[0].name);

            // Can we write empty strings in this range?
            REQUIRE(str_range[0] <= "");
            REQUIRE(str_range[1] >= "");
            // Can we write ASCII values in this range?
            REQUIRE(str_range[0] < " ");
            REQUIRE(str_range[1] > "~");

            std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(dim_infos[1].name);
            REQUIRE(u32_range[0] == (uint32_t)0);
            REQUIRE(u32_range[1] == (uint32_t)dim_infos[1].dim_max);

            // Check shape before write
            std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(!actual.has_value());

            // Check domainish accessors before resize
            ArrowTable non_empty_domain = sdf->get_non_empty_domain();
            std::vector<std::string> ned_str = ArrowAdapter::get_table_string_column_by_name(
                non_empty_domain, "mystring");

            ArrowTable soma_domain = sdf->get_soma_domain();
            std::vector<std::string> dom_str = ArrowAdapter::get_table_string_column_by_name(soma_domain, "mystring");

            ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
            std::vector<std::string> maxdom_str = ArrowAdapter::get_table_string_column_by_name(
                soma_maxdomain, "mystring");

            REQUIRE(ned_str == std::vector<std::string>({"", ""}));

            if (specify_domain) {
                REQUIRE(dom_str[0] == dim_infos[0].string_lo);
                REQUIRE(dom_str[1] == dim_infos[0].string_hi);
            } else {
                REQUIRE(dom_str == std::vector<std::string>({"", ""}));
            }
            REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));

            sdf->close();

            REQUIRE(sdf->nnz() == 0);

            // Write
            write_sjid_u32_str_data_from(0);

            REQUIRE(sdf->nnz() == 2);
            write_sjid_u32_str_data_from(8);
            // soma_joinid is not a dim here and so the second write is
            // an overwrite of the first here
            REQUIRE(sdf->nnz() == 2);

            // Check shape after write
            sdf = open(OpenMode::soma_read);
            actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(!actual.has_value());
            sdf->close();

            // Check can_upgrade_domain
            sdf = open(OpenMode::soma_read);
            StatusAndReason check = sdf->can_upgrade_soma_joinid_shape(1, "testing");
            // Must fail since this is too small.
            REQUIRE(check.first == false);
            REQUIRE(check.second == "testing: dataframe already has its domain set.");

            sdf->close();

            // Resize
            int64_t new_shape = int64_t{SOMA_JOINID_RESIZE_DIM_MAX + 1};
            uint32_t new_u32_dim_max = u32_dim_max * 2 + 1;

            // Check shape after write
            sdf = open(OpenMode::soma_read);
            actual = sdf->maybe_soma_joinid_shape();
            REQUIRE(!actual.has_value());
            sdf->close();

            // Apply the domain change
            if (test_upgrade_domain) {
                managed_unique_ptr<ArrowSchema> domain_schema = create_index_cols_info_schema(dim_infos);
                auto domain_array = ArrowAdapter::make_arrow_array_parent(dim_infos.size());
                domain_array->children[0] = ArrowAdapter::make_arrow_array_child_string(
                    std::vector<std::string>({"", ""}));
                domain_array->children[1] = ArrowAdapter::make_arrow_array_child(
                    std::vector<uint32_t>({0, new_u32_dim_max}));
                auto domain_table = ArrowTable(std::move(domain_array), std::move(domain_schema));

                // Not open for write
                sdf = open(OpenMode::soma_read);
                REQUIRE_THROWS(sdf->change_domain(domain_table, "testing"));
                sdf->close();

                // Open for write
                sdf = open(OpenMode::soma_write);
                sdf->change_domain(domain_table, "testing");
                sdf->close();
            } else {
                // Not open for write
                sdf = open(OpenMode::soma_read);
                REQUIRE_THROWS(sdf->resize_soma_joinid_shape(new_shape, "testing"));
                sdf->close();

                // Open for write
                sdf = open(OpenMode::soma_write);
                sdf->resize_soma_joinid_shape(new_shape, "testing");
                sdf->close();
            }

            sdf = open(OpenMode::soma_write);
            write_sjid_u32_str_data_from(SOMA_JOINID_DIM_MAX);
            sdf->close();

            // Check domainish accessors after resize
            sdf->open(OpenMode::soma_read, TimestampRange(0, 2));

            non_empty_domain = sdf->get_non_empty_domain();
            ned_str = ArrowAdapter::get_table_string_column_by_name(non_empty_domain, "mystring");

            soma_domain = sdf->get_soma_domain();
            dom_str = ArrowAdapter::get_table_string_column_by_name(soma_domain, "mystring");

            soma_maxdomain = sdf->get_soma_maxdomain();
            maxdom_str = ArrowAdapter::get_table_string_column_by_name(soma_maxdomain, "mystring");

            REQUIRE(ned_str == std::vector<std::string>({"", ""}));

            if (specify_domain) {
                REQUIRE(dom_str[0] == dim_infos[0].string_lo);
                REQUIRE(dom_str[1] == dim_infos[0].string_hi);
            } else {
                REQUIRE(dom_str == std::vector<std::string>({"", ""}));
            }
            REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));

            // Check can_resize_soma_joinid_shape
            check = sdf->can_resize_soma_joinid_shape(0, "testing");

            // Must pass since soma_joinid isn't a dim in this case.
            REQUIRE(check.first == true);
            REQUIRE(check.second == "");

            sdf->close();
        }
    }
}

TEST_CASE("SOMADataFrame: delete with only soma_joinid index", "[SOMADataFrame][delete]") {
    // Create a dataframe with soma_joinid as the only index column.
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://test-deletes-soma-joinid-only";

    {
        INFO("Create the dataframe.");
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(
            {helper::DimInfo(
                {.name = "soma_joinid",
                 .tiledb_datatype = TILEDB_INT64,
                 .dim_max = 7,
                 .string_lo = "N/A",
                 .string_hi = "N/A"})},
            {helper::AttrInfo({.name = "attr1", .tiledb_datatype = TILEDB_INT32})});
        SOMADataFrame::create(uri, schema, index_columns, ctx);
    }

    {
        INFO("Write data to array.");
        std::vector<int64_t> join(8);
        std::vector<int32_t> data(8);
        std::iota(join.begin(), join.end(), 0);
        std::iota(data.begin(), data.end(), 1);

        Array array{*ctx->tiledb_ctx(), uri, TILEDB_WRITE};
        Query query{*ctx->tiledb_ctx(), array};
        query.set_layout(TILEDB_GLOBAL_ORDER);
        query.set_data_buffer("soma_joinid", join);
        query.set_data_buffer("attr1", data);
        query.submit();
        query.finalize();
        REQUIRE(query.query_status() == tiledb::Query::Status::COMPLETE);
        array.close();
    }

    // Create variable for tests.
    int64_t expected_result_num{0};
    std::vector<int64_t> expected_joinids{};
    std::vector<int32_t> expected_data{};

    auto check_delete = [&](const std::string& log_note) {
        INFO(log_note);
        expected_joinids.resize(8, 0);
        expected_data.resize(8, 0);

        std::vector<int64_t> actual_joinids(8);
        std::vector<int32_t> actual_data(8);
        Array array{*ctx->tiledb_ctx(), uri, TILEDB_READ};
        Query query{*ctx->tiledb_ctx(), array};
        Subarray subarray(*ctx->tiledb_ctx(), array);
        subarray.add_range<int64_t>(0, 0, 7);
        query.set_layout(TILEDB_GLOBAL_ORDER);
        query.set_subarray(subarray);
        query.set_data_buffer("soma_joinid", actual_joinids);
        query.set_data_buffer("attr1", actual_data);
        query.submit();
        query.finalize();
        REQUIRE(query.query_status() == tiledb::Query::Status::COMPLETE);
        array.close();

        auto actual_result_num = static_cast<int64_t>(query.result_buffer_elements()["attr1"].second);
        CHECK(actual_result_num == expected_result_num);
        CHECK_THAT(actual_joinids, Catch::Matchers::Equals(expected_joinids));
        CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
    };

    SECTION("Error - must have some constraint") {
        INFO("Check throws if no constraint set.");
        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto delete_filter = df->create_coordinate_value_filter();
        CHECK_THROWS_AS(df->delete_cells(delete_filter), std::invalid_argument);
        df->close();
    }

    SECTION("Delete all by slice") {
        INFO("Delete all by slice");
        expected_result_num = 0;

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto delete_filter = df->create_coordinate_value_filter();
        delete_filter.add_slice<int64_t>(0, SliceSelection<int64_t>(-10, 10));
        df->delete_cells(delete_filter);
        df->close();

        check_delete("Delete all by slice");
    }

    SECTION("Delete all by points") {
        INFO("Delete all by points.");
        expected_result_num = 0;

        std::vector<int64_t> points{0, 3, 2, 1, 4, 7, 5, 6};
        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto delete_filter = df->create_coordinate_value_filter();
        delete_filter.add_points<int64_t>(0, PointSelection<int64_t>(points));
        df->delete_cells(delete_filter);
        df->close();

        check_delete("Delete all by points");
    }

    SECTION("Delete final value with slice") {
        expected_result_num = 7;
        expected_joinids.assign({0, 1, 2, 3, 4, 5, 6});
        expected_data.assign({1, 2, 3, 4, 5, 6, 7});

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto delete_filter = df->create_coordinate_value_filter();
        delete_filter.add_slice<int64_t>(0, SliceSelection<int64_t>(7, 10));
        df->delete_cells(delete_filter);
        df->close();

        check_delete("Delete final value with slice");
    }

    SECTION("Delete with coords and value filter") {
        expected_result_num = 6;
        expected_joinids.assign({0, 1, 2, 5, 6, 7});
        expected_data.assign({1, 2, 3, 6, 7, 8});

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_slice<int64_t>(0, SliceSelection<int64_t>(0, 4));
        int32_t max_value{3};
        auto value_filter = QueryCondition::create(*ctx->tiledb_ctx(), "attr1", max_value, TILEDB_GT);
        df->delete_cells(coord_filters, value_filter);
        df->close();

        check_delete("Delete using coords and value filter");
    }

    SECTION("Delete using value filter only") {
        expected_result_num = 4;
        expected_joinids.assign({0, 1, 2, 3});
        expected_data.assign({1, 2, 3, 4});

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        int32_t max_value{4};
        auto value_filter = QueryCondition::create(*ctx->tiledb_ctx(), "attr1", max_value, TILEDB_GT);
        df->delete_cells(coord_filters, value_filter);
        df->close();

        check_delete("Delete using value filter only");
    }
}

TEST_CASE("SOMADataFrame: delete with only string index column", "[SOMADataFrame][delete][string-index]") {
    // Create a dataframe with soma_joinid as the only index column.
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://test-deletes-string-index-only";

    {
        INFO("Create the dataframe.");
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(
            {helper::DimInfo(
                {.name = "label",
                 .tiledb_datatype = TILEDB_STRING_ASCII,
                 .dim_max = 0,
                 .string_lo = "",
                 .string_hi = ""})},
            {helper::AttrInfo({.name = "soma_joinid", .tiledb_datatype = TILEDB_INT64})});
        SOMADataFrame::create(uri, schema, index_columns, ctx);
    }

    std::vector<std::string> coords{"apple", "banana", "coconut", "durian", "eggplant", "fig"};
    uint64_t data_size = 0;
    for (auto& val : coords) {
        data_size += val.size();
    }
    std::vector<uint64_t> coords_offsets{};
    std::vector<char> coords_data(data_size);
    uint64_t curr_offset = 0;
    for (auto& val : coords) {
        coords_offsets.push_back(curr_offset);
        memcpy(coords_data.data() + curr_offset, val.data(), val.size());
        curr_offset += val.size();
    }

    std::vector<int64_t> index(6);
    std::iota(index.begin(), index.end(), 0);

    {
        INFO("Write data to array.");
        Array array{*ctx->tiledb_ctx(), uri, TILEDB_WRITE};
        Query query{*ctx->tiledb_ctx(), array};
        query.set_data_buffer("label", coords_data);
        query.set_offsets_buffer("label", coords_offsets);
        query.set_data_buffer("soma_joinid", index);
        query.submit();
        query.finalize();
        REQUIRE(query.query_status() == tiledb::Query::Status::COMPLETE);
    }

    // Create variable for tests.
    int64_t expected_result_num{0};
    std::vector<std::string> expected_labels{};
    std::vector<int64_t> expected_joinids{};

    auto check_delete = [&](const std::string& log_note) {
        INFO(log_note);

        // Create buffers for expected string labels.
        std::vector<char> expected_label_data(data_size);
        std::vector<uint64_t> expected_label_offsets{};
        uint64_t curr_offset = 0;
        for (auto& elem : expected_labels) {
            expected_label_offsets.push_back(curr_offset);
            memcpy(expected_label_data.data() + curr_offset, elem.data(), elem.size());
            curr_offset += elem.size();
        }
        expected_label_offsets.resize(7);
        expected_joinids.resize(6);

        std::vector<char> actual_label_data(data_size);
        std::vector<uint64_t> actual_label_offsets(7);
        std::vector<int64_t> actual_joinids(6);

        Array array{*ctx->tiledb_ctx(), uri, TILEDB_READ};
        Query query{*ctx->tiledb_ctx(), array};
        Subarray subarray(*ctx->tiledb_ctx(), array);
        subarray.add_range(0, std::string("a"), std::string("z"));
        query.set_subarray(subarray);
        query.set_layout(TILEDB_GLOBAL_ORDER);
        query.set_data_buffer("label", actual_label_data);
        query.set_offsets_buffer("label", actual_label_offsets);
        query.set_data_buffer("soma_joinid", actual_joinids);
        query.submit();
        query.finalize();
        REQUIRE(query.query_status() == tiledb::Query::Status::COMPLETE);
        REQUIRE(query.query_status() == tiledb::Query::Status::COMPLETE);
        array.close();

        auto actual_result_num = static_cast<int64_t>(query.result_buffer_elements()["soma_joinid"].second);
        CHECK(actual_result_num == expected_result_num);
        CHECK_THAT(actual_label_data, Catch::Matchers::Equals(expected_label_data));
        CHECK_THAT(actual_label_offsets, Catch::Matchers::Equals(expected_label_offsets));
        CHECK_THAT(actual_joinids, Catch::Matchers::Equals(expected_joinids));
    };

    SECTION("Delete all by slice") {
        expected_result_num = 0;

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_slice(0, SliceSelection<std::string>("a", "z"));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete all by slice");
    }

    SECTION("Delete all by points") {
        expected_result_num = 0;

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        std::vector<std::string> points{"fig", "banana", "coconut", "durian", "apple", "eggplant"};
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_points(0, PointSelection<std::string>(points));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete all by points");
    }

    SECTION("Delete with string slice") {
        expected_result_num = 3;
        expected_labels.assign({"apple", "eggplant", "fig"});
        expected_joinids.assign({0, 4, 5});

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_slice(0, SliceSelection<std::string>("b", "e"));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete with string slice");
    }

    SECTION("Delete with string points") {
        expected_result_num = 3;
        expected_labels.assign({"coconut", "eggplant", "fig"});
        expected_joinids.assign({2, 4, 5});
        std::vector<std::string> points{"banana", "durian", "apple"};

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_points(0, PointSelection<std::string>(points));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete with string points");
    }
}

TEST_CASE("SOMADataFrame: delete from multi-index dataframe", "[SOMADataFrame][delete]") {
    // Create a dataframe with soma_joinid as the only index column.
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://test-deletes-multi-index";

    {
        INFO("Create the dataframe.");
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(
            {helper::DimInfo(
                 {.name = "soma_joinid",
                  .tiledb_datatype = TILEDB_INT64,
                  .dim_max = 3,
                  .string_lo = "N/A",
                  .string_hi = "N/A"}),
             helper::DimInfo(
                 {.name = "index",
                  .tiledb_datatype = TILEDB_UINT32,
                  .dim_max = 2,
                  .string_lo = "N/A",
                  .string_hi = "N/A"})},
            {helper::AttrInfo({.name = "attr1", .tiledb_datatype = TILEDB_INT32})});
        SOMADataFrame::create(uri, schema, index_columns, ctx);
    }

    {
        INFO("Write data to array.");
        std::vector<int64_t> join{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
        std::vector<uint32_t> index{0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
        std::vector<int32_t> data(12);
        std::iota(data.begin(), data.end(), 1);

        Array array{*ctx->tiledb_ctx(), uri, TILEDB_WRITE};
        Query query{*ctx->tiledb_ctx(), array};
        query.set_data_buffer("soma_joinid", join);
        query.set_data_buffer("index", index);
        query.set_data_buffer("attr1", data);
        query.submit();
        query.finalize();
        REQUIRE(query.query_status() == tiledb::Query::Status::COMPLETE);
        array.close();
    }

    // Create variable for tests. Using sections will rerun the this test from beginning to end for each section.
    int64_t expected_result_num{0};
    std::vector<int64_t> expected_joinids{};
    std::vector<uint32_t> expected_index{};
    std::vector<int32_t> expected_data{};

    auto check_delete = [&](const std::string& log_note) {
        INFO(log_note);
        expected_joinids.resize(12, 0);
        expected_index.resize(12, 0);
        expected_data.resize(12, 0);

        std::vector<int64_t> actual_joinids(12);
        std::vector<uint32_t> actual_index(12);
        std::vector<int32_t> actual_data(12);
        Array array{*ctx->tiledb_ctx(), uri, TILEDB_READ};
        Query query{*ctx->tiledb_ctx(), array};
        Subarray subarray(*ctx->tiledb_ctx(), array);
        subarray.add_range<int64_t>(0, 0, 3);
        subarray.add_range<uint32_t>(1, 0, 2);
        query.set_subarray(subarray);
        query.set_data_buffer("soma_joinid", actual_joinids);
        query.set_data_buffer("index", actual_index);
        query.set_data_buffer("attr1", actual_data);
        query.submit();
        query.finalize();
        REQUIRE(query.query_status() == tiledb::Query::Status::COMPLETE);
        array.close();

        auto actual_result_num = static_cast<int64_t>(query.result_buffer_elements()["attr1"].second);
        CHECK(actual_result_num == expected_result_num);
        CHECK_THAT(actual_joinids, Catch::Matchers::Equals(expected_joinids));
        CHECK_THAT(actual_index, Catch::Matchers::Equals(expected_index));
        CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
    };

    SECTION("Delete all by joinid slice") {
        expected_result_num = 0;

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_slice<int64_t>(0, SliceSelection<int64_t>(0, 3));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete all by joinid slice");
    }

    SECTION("Delete all by joinid slice and monostate") {
        expected_result_num = 0;

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_column_selection<int64_t>(0, SliceSelection<int64_t>(0, 3));
        coord_filters.add_column_selection<uint32_t>(1, std::monostate());
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete all by joinid slice and monostate");
    }

    SECTION("Delete all by joinid points") {
        expected_result_num = 0;

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        std::vector<int64_t> points{1, 3, 2, 0};
        coord_filters.add_points<int64_t>(0, PointSelection<int64_t>(points));
        df->delete_cells(coord_filters);
        df->close();
        check_delete("Delete all by joinid points");
    }

    SECTION("Delete all by index slice") {
        expected_result_num = 0;

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_slice<uint32_t>(1, SliceSelection<uint32_t>(0, 2));
        df->delete_cells(coord_filters);
        df->close();
        check_delete("Delete all by joinid slice");
    }

    SECTION("Delete all by index points") {
        expected_result_num = 0;

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        std::vector<uint32_t> points{2, 0, 1};
        coord_filters.add_points<uint32_t>(1, PointSelection<uint32_t>(points));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete all by joinid slice");
    }
    SECTION("Delete by joinid slice") {
        expected_result_num = 6;
        expected_joinids.assign({0, 0, 0, 1, 1, 1});
        expected_index.assign({0, 1, 2, 0, 1, 2});
        expected_data.assign({1, 2, 3, 4, 5, 6});

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_slice<int64_t>(0, SliceSelection<int64_t>(2, 10));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete by joinid slice");
    }
    SECTION("Delete by joinid slice and all indices") {
        expected_result_num = 6;
        expected_joinids.assign({0, 0, 0, 1, 1, 1});
        expected_index.assign({0, 1, 2, 0, 1, 2});
        expected_data.assign({1, 2, 3, 4, 5, 6});

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_slice(0, SliceSelection<int64_t>(2, 10));
        coord_filters.add_slice(1, SliceSelection<uint32_t>(0, 10));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete by joinid slice");
    }
    SECTION("Delete multiple slice selections") {
        expected_result_num = 8;
        expected_joinids.assign({0, 0, 0, 1, 2, 3, 3, 3});
        expected_index.assign({0, 1, 2, 2, 2, 0, 1, 2});
        expected_data.assign({1, 2, 3, 6, 9, 10, 11, 12});

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_slice(1, SliceSelection<uint32_t>(0, 1));
        coord_filters.add_slice(0, SliceSelection<int64_t>(1, 2));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete multiple slice selections");
    }
    SECTION("Delete multiple point selections") {
        expected_result_num = 10;
        expected_joinids.assign({0, 0, 0, 1, 1, 2, 2, 3, 3, 3});
        expected_index.assign({0, 1, 2, 0, 1, 0, 1, 0, 1, 2});
        expected_data.assign({1, 2, 3, 4, 5, 7, 8, 10, 11, 12});
        std::vector<int64_t> join_points{2, 1};
        std::vector<uint32_t> index_points{2};

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_points(1, PointSelection<uint32_t>(index_points));
        coord_filters.add_points(0, PointSelection<int64_t>(join_points));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete multiple point selections");
    }
    SECTION("Delete mixed selection") {
        expected_result_num = 10;
        expected_joinids.assign({0, 0, 0, 1, 1, 2, 2, 3, 3, 3});
        expected_index.assign({0, 1, 2, 0, 1, 0, 1, 0, 1, 2});
        expected_data.assign({1, 2, 3, 4, 5, 7, 8, 10, 11, 12});
        std::vector<uint32_t> index_points{2};

        auto df = SOMADataFrame::open(uri, OpenMode::soma_delete, ctx, std::nullopt);
        auto coord_filters = df->create_coordinate_value_filter();
        coord_filters.add_points(1, PointSelection<uint32_t>(index_points));
        coord_filters.add_slice(0, SliceSelection<int64_t>(1, 2));
        df->delete_cells(coord_filters);
        df->close();

        check_delete("Delete mixed selections");
    }
}
