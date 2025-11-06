/**
 * @file   unit_soma_column.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for implementation of SOMAColumn class
 */

#include <tiledb/tiledb>
#include <tiledbsoma/tiledbsoma>
#include <tuple>
#include "../src/soma/soma_attribute.h"
#include "../src/soma/soma_column.h"
#include "../src/soma/soma_dimension.h"
#include "../src/soma/soma_geometry_column.h"
#include "common.h"

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
    static const inline int64_t i64_dim_max = 99;
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

TEST_CASE("SOMAColumn: SOMADimension") {
    auto ctx = std::make_shared<SOMAContext>();
    PlatformConfig platform_config{};

    std::vector<helper::DimInfo> dim_infos(
        {helper::DimInfo(
             {.name = "dimension",
              .tiledb_datatype = TILEDB_UINT32,
              .dim_max = 100,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "dimension",
              .tiledb_datatype = TILEDB_FLOAT64,
              .dim_max = 100,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "dimension",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = 100,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "dimension",
              .tiledb_datatype = TILEDB_STRING_ASCII,
              .dim_max = 100,
              .string_lo = "N/A",
              .string_hi = "N/A"})});

    std::vector<helper::DimInfo> geom_dim_infos({helper::DimInfo(
        {.name = "dimension",
         .tiledb_datatype = TILEDB_GEOM_WKB,
         .dim_max = 100,
         .string_lo = "N/A",
         .string_hi = "N/A"})});

    std::vector<helper::DimInfo> spatial_dim_infos(
        {helper::DimInfo(
             {.name = "x", .tiledb_datatype = TILEDB_FLOAT64, .dim_max = 200, .string_lo = "N/A", .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "y",
              .tiledb_datatype = TILEDB_FLOAT64,
              .dim_max = 100,
              .string_lo = "N/A",
              .string_hi = "N/A"})});

    auto index_columns = helper::create_column_index_info(dim_infos);

    std::vector<std::shared_ptr<SOMAColumn>> columns;

    for (int64_t i = 0; i < index_columns.second->n_children; ++i) {
        columns.push_back(
            SOMADimension::create(
                ctx->tiledb_ctx(),
                index_columns.second->children[i],
                index_columns.first->children[i],
                "SOMAGeometryDataFrame",
                "",
                platform_config));

        REQUIRE(columns.back()->tiledb_dimensions().value()[0].type() == dim_infos[i].tiledb_datatype);
    }

    REQUIRE(columns[1]->core_domain_slot<double_t>() == std::make_pair<double_t, double_t>(0, helper::CORE_DOMAIN_MAX));
    REQUIRE(columns[1]->core_domain_slot<double_t>() == std::make_pair<double_t, double_t>(0, helper::CORE_DOMAIN_MAX));
    REQUIRE(columns[2]->core_domain_slot<int64_t>() == std::make_pair<int64_t, int64_t>(0, helper::CORE_DOMAIN_MAX));
    REQUIRE(columns[3]->core_domain_slot<std::string>() == std::make_pair<std::string, std::string>("", ""));
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMAColumn: query variant-indexed dataframe dim-str-u32 attr-sjid",
    "[SOMAColumn]") {
    auto specify_domain = GENERATE(false, true);
    SECTION("- specify_domain=" + std::to_string(specify_domain)) {
        std::string suffix1 = specify_domain ? "true" : "false";
        set_up(std::make_shared<SOMAContext>(), "mem://unit-test-column-variant-indexed-dataframe-4-" + suffix1);

        std::string string_lo = "";
        std::string string_hi = "";
        std::vector<helper::DimInfo> dim_infos({str_dim_info(string_lo, string_hi), u32_dim_info()});
        std::vector<helper::AttrInfo> attr_infos({i64_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::soma_read);

        // External column initialization
        auto raw_array = tiledb::Array(*ctx_->tiledb_ctx(), uri_, TILEDB_READ);
        std::vector<std::shared_ptr<SOMAColumn>> columns;

        for (auto dimension : sdf->tiledb_schema()->domain().dimensions()) {
            columns.push_back(std::make_shared<SOMADimension>(SOMADimension(dimension)));
        }

        for (size_t i = 0; i < sdf->tiledb_schema()->attribute_num(); ++i) {
            columns.push_back(std::make_shared<SOMAAttribute>(SOMAAttribute(sdf->tiledb_schema()->attribute(i))));
        }

        CurrentDomain current_domain = sdf->get_current_domain_for_test();

        REQUIRE(!current_domain.is_empty());
        REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
        NDRectangle ndrect = current_domain.ndrectangle();

        std::array<std::string, 2> str_range = ndrect.range<std::string>(dim_infos[0].name);
        std::pair<std::string, std::string> str_external = columns[0]->core_current_domain_slot<std::string>(
            *ctx_, raw_array);

        // Can we write empty strings in this range?
        REQUIRE(str_range[0] <= "");
        REQUIRE(str_range[1] >= "");
        // Can we write ASCII values in this range?
        REQUIRE(str_range[0] < " ");
        REQUIRE(str_range[1] > "~");

        // SOMAColumn return the current domain to return to the user directly
        // ("", "") and not the internal values (e.g. \xff or \x7f)
        REQUIRE(str_external.first == "");
        REQUIRE(str_external.second == "");

        std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(dim_infos[1].name);
        std::pair<uint32_t, uint32_t> u32_external = columns[1]->core_current_domain_slot<uint32_t>(*ctx_, raw_array);
        REQUIRE(u32_range[0] == u32_external.first);
        REQUIRE(u32_range[1] == u32_external.second);

        // Check shape before write
        std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(!actual.has_value());

        // Check domainish accessors before resize
        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<std::string> ned_str = ArrowAdapter::get_table_string_column_by_name(non_empty_domain, "mystring");

        auto col_non_empty_domain = columns[0]->arrow_domain_slot(*ctx_, raw_array, Domainish::kind_non_empty_domain);
        std::vector<std::string> ned_str_col = std::apply(ArrowAdapter::get_array_string_column, col_non_empty_domain);

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<std::string> dom_str = ArrowAdapter::get_table_string_column_by_name(soma_domain, "mystring");
        auto col_soma_domain = columns[0]->arrow_domain_slot(*ctx_, raw_array, Domainish::kind_core_current_domain);
        std::vector<std::string> dom_str_col = std::apply(ArrowAdapter::get_array_string_column, col_soma_domain);

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<std::string> maxdom_str = ArrowAdapter::get_table_string_column_by_name(soma_maxdomain, "mystring");

        auto col_soma_maxdomain = columns[0]->arrow_domain_slot(*ctx_, raw_array, Domainish::kind_core_domain);
        std::vector<std::string> maxdom_str_col = std::apply(ArrowAdapter::get_array_string_column, col_soma_maxdomain);

        // Cleanup domain arrow tables
        col_non_empty_domain.first->release(col_non_empty_domain.first);
        col_non_empty_domain.second->release(col_non_empty_domain.second);
        col_soma_domain.first->release(col_soma_domain.first);
        col_soma_domain.second->release(col_soma_domain.second);
        col_soma_maxdomain.first->release(col_soma_maxdomain.first);
        col_soma_maxdomain.second->release(col_soma_maxdomain.second);

        free(col_non_empty_domain.first);
        free(col_non_empty_domain.second);
        free(col_soma_domain.first);
        free(col_soma_domain.second);
        free(col_soma_maxdomain.first);
        free(col_soma_maxdomain.second);

        REQUIRE(ned_str == std::vector<std::string>({"", ""}));

        REQUIRE(ned_str == ned_str_col);
        REQUIRE(dom_str == dom_str_col);
        REQUIRE(maxdom_str == maxdom_str_col);

        if (specify_domain) {
            REQUIRE(dom_str[0] == dim_infos[0].string_lo);
            REQUIRE(dom_str[1] == dim_infos[0].string_hi);
        } else {
            REQUIRE(dom_str == std::vector<std::string>({"", ""}));
        }
        REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));

        sdf->close();

        sdf = open(OpenMode::soma_write);
        write_sjid_u32_str_data_from(0);

        sdf->close();

        sdf = open(OpenMode::soma_read);
        REQUIRE(sdf->nnz() == 2);

        sdf->close();

        auto external_query = open(OpenMode::soma_read)->create_managed_query();

        columns[1]->select_columns(external_query);
        columns[1]->set_dim_point<uint32_t>(external_query, 1234);

        // Configure query and allocate result buffers
        auto ext_res = external_query.read_next().value();

        REQUIRE(ext_res->num_rows() == 1);

        external_query.reset();

        const std::vector a({std::make_pair<std::string, std::string>("apple", "b")});

        columns[0]->select_columns(external_query);
        columns[0]->set_dim_ranges(
            external_query, std::vector({std::make_pair<std::string, std::string>("apple", "b")}));

        // Configure query and allocate result buffers
        ext_res = external_query.read_next().value();

        REQUIRE(ext_res->num_rows() == 1);
    }
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMAColumn: query variant-indexed dataframe dim-str-u32 attr-sjid",
    "[SOMADataFrame]") {
    auto specify_domain = GENERATE(false, true);
    SECTION("- specify_domain=" + std::to_string(specify_domain)) {
        std::string suffix1 = specify_domain ? "true" : "false";
        set_up(std::make_shared<SOMAContext>(), "mem://unit-test-column-variant-indexed-dataframe-4-" + suffix1);

        std::string string_lo = "";
        std::string string_hi = "";
        std::vector<helper::DimInfo> dim_infos({str_dim_info(string_lo, string_hi), u32_dim_info()});
        std::vector<helper::AttrInfo> attr_infos({i64_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::soma_read);

        // External column initialization
        auto raw_array = tiledb::Array(*ctx_->tiledb_ctx(), uri_, TILEDB_READ);
        std::vector<std::shared_ptr<SOMAColumn>> columns;

        for (auto dimension : sdf->tiledb_schema()->domain().dimensions()) {
            columns.push_back(std::make_shared<SOMADimension>(SOMADimension(dimension)));
        }

        for (size_t i = 0; i < sdf->tiledb_schema()->attribute_num(); ++i) {
            columns.push_back(std::make_shared<SOMAAttribute>(SOMAAttribute(sdf->tiledb_schema()->attribute(i))));
        }

        CurrentDomain current_domain = sdf->get_current_domain_for_test();

        REQUIRE(!current_domain.is_empty());
        REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
        NDRectangle ndrect = current_domain.ndrectangle();

        std::array<std::string, 2> str_range = ndrect.range<std::string>(dim_infos[0].name);
        std::pair<std::string, std::string> str_external = columns[0]->core_current_domain_slot<std::string>(
            *ctx_, raw_array);

        // Can we write empty strings in this range?
        REQUIRE(str_range[0] <= "");
        REQUIRE(str_external.first <= "");
        REQUIRE(str_range[1] >= "");
        REQUIRE(str_external.second >= "");
        // Can we write ASCII values in this range?
        REQUIRE(str_range[0] < " ");
        REQUIRE(str_external.first <= " ");
        REQUIRE(str_range[1] > "~");
        // REQUIRE(str_external.second >= "~");

        std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(dim_infos[1].name);
        std::pair<uint32_t, uint32_t> u32_external = columns[1]->core_current_domain_slot<uint32_t>(*ctx_, raw_array);
        REQUIRE(u32_range[0] == u32_external.first);
        REQUIRE(u32_range[1] == u32_external.second);

        // Check shape before write
        std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(!actual.has_value());

        // Check domainish accessors before resize
        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<std::string> ned_str = ArrowAdapter::get_table_string_column_by_name(non_empty_domain, "mystring");

        auto col_non_empty_domain = columns[0]->arrow_domain_slot(*ctx_, raw_array, Domainish::kind_non_empty_domain);
        std::vector<std::string> ned_str_col = std::apply(ArrowAdapter::get_array_string_column, col_non_empty_domain);

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<std::string> dom_str = ArrowAdapter::get_table_string_column_by_name(soma_domain, "mystring");

        auto col_soma_domain = columns[0]->arrow_domain_slot(*ctx_, raw_array, Domainish::kind_core_current_domain);
        std::vector<std::string> dom_str_col = std::apply(ArrowAdapter::get_array_string_column, col_soma_domain);

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<std::string> maxdom_str = ArrowAdapter::get_table_string_column_by_name(soma_maxdomain, "mystring");

        auto col_soma_maxdomain = columns[0]->arrow_domain_slot(*ctx_, raw_array, Domainish::kind_core_domain);
        std::vector<std::string> maxdom_str_col = std::apply(ArrowAdapter::get_array_string_column, col_soma_maxdomain);

        // Cleanup domain arrow tables
        col_non_empty_domain.first->release(col_non_empty_domain.first);
        col_non_empty_domain.second->release(col_non_empty_domain.second);
        col_soma_domain.first->release(col_soma_domain.first);
        col_soma_domain.second->release(col_soma_domain.second);
        col_soma_maxdomain.first->release(col_soma_maxdomain.first);
        col_soma_maxdomain.second->release(col_soma_maxdomain.second);

        free(col_non_empty_domain.first);
        free(col_non_empty_domain.second);
        free(col_soma_domain.first);
        free(col_soma_domain.second);
        free(col_soma_maxdomain.first);
        free(col_soma_maxdomain.second);

        REQUIRE(ned_str == std::vector<std::string>({"", ""}));

        REQUIRE(ned_str == ned_str_col);
        REQUIRE(dom_str == dom_str_col);
        REQUIRE(maxdom_str == maxdom_str_col);

        if (specify_domain) {
            REQUIRE(dom_str[0] == dim_infos[0].string_lo);
            REQUIRE(dom_str[1] == dim_infos[0].string_hi);
        } else {
            REQUIRE(dom_str == std::vector<std::string>({"", ""}));
        }
        REQUIRE(maxdom_str == std::vector<std::string>({"", ""}));

        sdf->close();

        sdf = open(OpenMode::soma_write);
        write_sjid_u32_str_data_from(0);

        sdf->close();

        sdf = open(OpenMode::soma_read);
        REQUIRE(sdf->nnz() == 2);

        sdf->close();

        auto external_query = open(OpenMode::soma_read)->create_managed_query();

        columns[1]->select_columns(external_query);
        columns[1]->set_dim_point<uint32_t>(external_query, 1234);

        // Configure query and allocate result buffers
        auto ext_res = external_query.read_next().value();

        REQUIRE(ext_res->num_rows() == 1);

        external_query.reset();

        columns[0]->select_columns(external_query);
        columns[0]->set_dim_ranges<std::string>(
            external_query, std::vector({std::make_pair<std::string, std::string>("apple", "b")}));

        // Configure query and allocate result buffers
        ext_res = external_query.read_next().value();

        REQUIRE(ext_res->num_rows() == 1);
    }
}

// Test enumeration handling in SOMAColumn::deserialize
TEST_CASE("SOMAColumn: Enumeration handling in deserialize") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-enumeration-handling";

    // Create array with enumeration
    auto vfs = VFS(*ctx->tiledb_ctx());
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    // Create schema with enumerated attribute
    ArraySchema schema(*ctx->tiledb_ctx(), TILEDB_SPARSE);
    auto dim = Dimension::create<int64_t>(*ctx->tiledb_ctx(), "d0", {0, 1000}, 10);
    Domain domain(*ctx->tiledb_ctx());
    domain.add_dimension(dim);
    schema.set_domain(domain);

    // Create enumeration for categorical data
    std::vector<std::string> enum_values = {"red", "green", "blue"};
    auto enumeration = Enumeration::create(*ctx->tiledb_ctx(), "color_enum", enum_values);
    ArraySchemaExperimental::add_enumeration(*ctx->tiledb_ctx(), schema, enumeration);

    // Create attribute with enumeration
    auto attr = Attribute::create<uint8_t>(*ctx->tiledb_ctx(), "color");
    AttributeExperimental::set_enumeration_name(*ctx->tiledb_ctx(), attr, "color_enum");
    schema.add_attribute(attr);

    // Create another attribute without enumeration
    auto attr2 = Attribute::create<int32_t>(*ctx->tiledb_ctx(), "value");
    schema.add_attribute(attr2);

    schema.check();
    Array::create(uri, schema);

    // Test deserialize with enumeration
    {
        Array array(*ctx->tiledb_ctx(), uri, TILEDB_READ);
        auto columns = SOMAColumn::deserialize(*ctx->tiledb_ctx(), array, {});

        // Should have 3 columns: 1 dimension + 2 attributes
        REQUIRE(columns.size() == 3);

        // Check that enumerated attribute is properly handled
        bool found_enumerated = false;
        bool found_non_enumerated = false;

        for (const auto& column : columns) {
            if (column->tiledb_attributes().has_value()) {
                auto attributes = column->tiledb_attributes().value();
                for (const auto& attribute : attributes) {
                    if (attribute.name() == "color") {
                        // Check that enumeration is detected
                        auto enum_name = AttributeExperimental::get_enumeration_name(*ctx->tiledb_ctx(), attribute);
                        REQUIRE(enum_name.has_value());
                        REQUIRE(enum_name.value() == "color_enum");
                        found_enumerated = true;
                    } else if (attribute.name() == "value") {
                        // Check that no enumeration for regular attribute
                        auto enum_name = AttributeExperimental::get_enumeration_name(*ctx->tiledb_ctx(), attribute);
                        REQUIRE(!enum_name.has_value());
                        found_non_enumerated = true;
                    }
                }
            }
        }

        REQUIRE(found_enumerated);
        REQUIRE(found_non_enumerated);
        array.close();
    }
}

// Test enumeration handling when metadata is missing
TEST_CASE("SOMAColumn: Enumeration handling without metadata") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-enumeration-no-metadata";

    // Create array with enumeration but no SOMA metadata
    auto vfs = VFS(*ctx->tiledb_ctx());
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    ArraySchema schema(*ctx->tiledb_ctx(), TILEDB_SPARSE);
    auto dim = Dimension::create<int64_t>(*ctx->tiledb_ctx(), "d0", {0, 1000}, 10);
    Domain domain(*ctx->tiledb_ctx());
    domain.add_dimension(dim);
    schema.set_domain(domain);

    // Create enumeration
    std::vector<std::string> enum_values = {"cat", "dog", "bird"};
    auto enumeration = Enumeration::create(*ctx->tiledb_ctx(), "animal_enum", enum_values);
    ArraySchemaExperimental::add_enumeration(*ctx->tiledb_ctx(), schema, enumeration);

    // Create attribute with enumeration
    auto attr = Attribute::create<uint8_t>(*ctx->tiledb_ctx(), "animal");
    AttributeExperimental::set_enumeration_name(*ctx->tiledb_ctx(), attr, "animal_enum");
    schema.add_attribute(attr);

    schema.check();
    Array::create(uri, schema);

    // Test deserialize without SOMA metadata (empty metadata map)
    {
        Array array(*ctx->tiledb_ctx(), uri, TILEDB_READ);
        std::map<std::string, tiledbsoma::MetadataValue> empty_metadata;
        auto columns = SOMAColumn::deserialize(*ctx->tiledb_ctx(), array, empty_metadata);

        // Should have 2 columns: 1 dimension + 1 attribute
        REQUIRE(columns.size() == 2);

        // Check that enumerated attribute is still properly handled
        bool found_enumerated = false;

        for (const auto& column : columns) {
            if (column->tiledb_attributes().has_value()) {
                auto attributes = column->tiledb_attributes().value();
                for (const auto& attribute : attributes) {
                    if (attribute.name() == "animal") {
                        auto enum_name = AttributeExperimental::get_enumeration_name(*ctx->tiledb_ctx(), attribute);
                        REQUIRE(enum_name.has_value());
                        REQUIRE(enum_name.value() == "animal_enum");
                        found_enumerated = true;
                    }
                }
            }
        }

        REQUIRE(found_enumerated);
        array.close();
    }
}

// Test error cases in core_current_domain_slot for strings
TEST_CASE("SOMAColumn: String domain error handling") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-string-domain-errors";

    // Create array with string dimension
    auto vfs = VFS(*ctx->tiledb_ctx());
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    ArraySchema schema(*ctx->tiledb_ctx(), TILEDB_SPARSE);
    auto dim = Dimension::create(*ctx->tiledb_ctx(), "str_dim", TILEDB_STRING_ASCII, nullptr, nullptr);
    dim.set_cell_val_num(TILEDB_VAR_NUM);
    Domain domain(*ctx->tiledb_ctx());
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int32_t>(*ctx->tiledb_ctx(), "value");
    schema.add_attribute(attr);
    schema.check();

    Array::create(uri, schema);

    // Test core_current_domain_slot with valid string ranges
    {
        Array array(*ctx->tiledb_ctx(), uri, TILEDB_READ);
        auto columns = SOMAColumn::deserialize(*ctx->tiledb_ctx(), array, {});

        // Find the string dimension column
        std::shared_ptr<SOMAColumn> string_dim_column = nullptr;
        for (const auto& column : columns) {
            if (column->tiledb_dimensions().has_value()) {
                auto dimensions = column->tiledb_dimensions().value();
                for (const auto& dimension : dimensions) {
                    if (dimension.name() == "str_dim") {
                        string_dim_column = column;
                        break;
                    }
                }
            }
        }

        REQUIRE(string_dim_column != nullptr);

        // Test core_domain_slot for strings (should return empty strings)
        auto core_domain = string_dim_column->core_domain_slot<std::string>();
        REQUIRE(core_domain.first == "");
        REQUIRE(core_domain.second == "");

        array.close();
    }
}

// Test current domain slot error handling with SOMADataFrame
TEST_CASE("SOMAColumn: String current domain with SOMADataFrame") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-string-current-domain-sdf";

    // Create SOMADataFrame with string dimension (this automatically sets up current domain)
    std::vector<helper::DimInfo> dim_infos({helper::DimInfo(
        {.name = "str_dim",
         .tiledb_datatype = TILEDB_STRING_ASCII,
         .dim_max = 0,  // Not used for string dims
         .string_lo = "",
         .string_hi = ""})});

    std::vector<helper::AttrInfo> attr_infos({helper::AttrInfo({.name = "value", .tiledb_datatype = TILEDB_INT32})});

    auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(dim_infos, attr_infos);
    SOMADataFrame::create(uri, schema, index_columns, ctx);

    // Test current domain slot using SOMADataFrame approach
    {
        auto sdf = SOMADataFrame::open(uri, OpenMode::soma_read, ctx);
        auto raw_array = tiledb::Array(*ctx->tiledb_ctx(), uri, TILEDB_READ);
        auto columns = SOMAColumn::deserialize(*ctx->tiledb_ctx(), raw_array, {});

        // Find the string dimension column
        std::shared_ptr<SOMAColumn> string_dim_column = nullptr;
        for (const auto& column : columns) {
            if (column->tiledb_dimensions().has_value()) {
                auto dimensions = column->tiledb_dimensions().value();
                for (const auto& dimension : dimensions) {
                    if (dimension.name() == "str_dim") {
                        string_dim_column = column;
                        break;
                    }
                }
            }
        }

        REQUIRE(string_dim_column != nullptr);

        // Test current domain functionality with SOMADataFrame
        if (sdf->has_current_domain()) {
            auto current_domain_pair = string_dim_column->core_current_domain_slot<std::string>(*ctx, raw_array);
            REQUIRE(current_domain_pair.first == "");
            REQUIRE(current_domain_pair.second == "");
        }

        sdf->close();
        raw_array.close();
    }
}

// Test current domain slot error path coverage
TEST_CASE("SOMAColumn: String current domain error path coverage") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-string-current-domain-errors";

    // Create SOMADataFrame with string dimension
    std::vector<helper::DimInfo> dim_infos({helper::DimInfo(
        {.name = "str_dim", .tiledb_datatype = TILEDB_STRING_ASCII, .dim_max = 0, .string_lo = "", .string_hi = ""})});

    std::vector<helper::AttrInfo> attr_infos({helper::AttrInfo({.name = "value", .tiledb_datatype = TILEDB_INT32})});

    auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(dim_infos, attr_infos);
    SOMADataFrame::create(uri, schema, index_columns, ctx);

    // Test error handling paths in current domain methods
    {
        auto sdf = SOMADataFrame::open(uri, OpenMode::soma_read, ctx);
        auto raw_array = tiledb::Array(*ctx->tiledb_ctx(), uri, TILEDB_READ);
        auto columns = SOMAColumn::deserialize(*ctx->tiledb_ctx(), raw_array, {});

        // Find the string dimension column
        std::shared_ptr<SOMAColumn> string_dim_column = nullptr;
        for (const auto& column : columns) {
            if (column->tiledb_dimensions().has_value()) {
                auto dimensions = column->tiledb_dimensions().value();
                for (const auto& dimension : dimensions) {
                    if (dimension.name() == "str_dim") {
                        string_dim_column = column;
                        break;
                    }
                }
            }
        }

        REQUIRE(string_dim_column != nullptr);

        // Test basic current domain functionality
        auto current_domain_pair = string_dim_column->core_current_domain_slot<std::string>(*ctx, raw_array);
        REQUIRE(current_domain_pair.first == "");
        REQUIRE(current_domain_pair.second == "");

        // Test with valid current domain if available
        if (sdf->has_current_domain()) {
            auto current_domain = sdf->get_current_domain_for_test();
            if (!current_domain.is_empty() && current_domain.type() == TILEDB_NDRECTANGLE) {
                NDRectangle ndrect = current_domain.ndrectangle();
                auto ndrect_domain_pair = string_dim_column->core_current_domain_slot<std::string>(ndrect);
                REQUIRE(ndrect_domain_pair.first == "");
                REQUIRE(ndrect_domain_pair.second == "");
            }
        }

        sdf->close();
        raw_array.close();
    }
}

// Test multiple enumerated attributes in a single array
TEST_CASE("SOMAColumn: Multiple enumerated attributes") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-multiple-enumerations";

    // Create array with multiple enumerated attributes
    auto vfs = VFS(*ctx->tiledb_ctx());
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    ArraySchema schema(*ctx->tiledb_ctx(), TILEDB_SPARSE);
    auto dim = Dimension::create<int64_t>(*ctx->tiledb_ctx(), "d0", {0, 1000}, 10);
    Domain domain(*ctx->tiledb_ctx());
    domain.add_dimension(dim);
    schema.set_domain(domain);

    // Create first enumeration
    std::vector<std::string> color_values = {"red", "green", "blue"};
    auto color_enum = Enumeration::create(*ctx->tiledb_ctx(), "color_enum", color_values);
    ArraySchemaExperimental::add_enumeration(*ctx->tiledb_ctx(), schema, color_enum);

    // Create second enumeration
    std::vector<std::string> size_values = {"small", "medium", "large"};
    auto size_enum = Enumeration::create(*ctx->tiledb_ctx(), "size_enum", size_values);
    ArraySchemaExperimental::add_enumeration(*ctx->tiledb_ctx(), schema, size_enum);

    // Create attributes with different enumerations
    auto attr1 = Attribute::create<uint8_t>(*ctx->tiledb_ctx(), "color");
    AttributeExperimental::set_enumeration_name(*ctx->tiledb_ctx(), attr1, "color_enum");
    schema.add_attribute(attr1);

    auto attr2 = Attribute::create<uint8_t>(*ctx->tiledb_ctx(), "size");
    AttributeExperimental::set_enumeration_name(*ctx->tiledb_ctx(), attr2, "size_enum");
    schema.add_attribute(attr2);

    // Add a non-enumerated attribute
    auto attr3 = Attribute::create<int32_t>(*ctx->tiledb_ctx(), "value");
    schema.add_attribute(attr3);

    schema.check();
    Array::create(uri, schema);

    // Test deserialization with multiple enumerations
    {
        Array array(*ctx->tiledb_ctx(), uri, TILEDB_READ);
        auto columns = SOMAColumn::deserialize(*ctx->tiledb_ctx(), array, {});

        // Should have 4 columns: 1 dimension + 3 attributes
        REQUIRE(columns.size() == 4);

        std::map<std::string, std::string> found_enumerations;
        bool found_non_enumerated = false;

        for (const auto& column : columns) {
            if (column->tiledb_attributes().has_value()) {
                auto attributes = column->tiledb_attributes().value();
                for (const auto& attribute : attributes) {
                    auto enum_name = AttributeExperimental::get_enumeration_name(*ctx->tiledb_ctx(), attribute);
                    if (enum_name.has_value()) {
                        found_enumerations[attribute.name()] = enum_name.value();
                    } else if (attribute.name() == "value") {
                        found_non_enumerated = true;
                    }
                }
            }
        }

        REQUIRE(found_enumerations.size() == 2);
        REQUIRE(found_enumerations["color"] == "color_enum");
        REQUIRE(found_enumerations["size"] == "size_enum");
        REQUIRE(found_non_enumerated);

        array.close();
    }
}

// Test NDRectangle version using existing working pattern
TEST_CASE("SOMAColumn: String current domain NDRectangle functionality") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-string-ndrectangle-func";

    // Create SOMADataFrame with string dimension using working pattern
    std::vector<helper::DimInfo> dim_infos({helper::DimInfo(
        {.name = "str_dim", .tiledb_datatype = TILEDB_STRING_ASCII, .dim_max = 0, .string_lo = "", .string_hi = ""})});

    std::vector<helper::AttrInfo> attr_infos({helper::AttrInfo({.name = "value", .tiledb_datatype = TILEDB_INT32})});

    auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(dim_infos, attr_infos);
    SOMADataFrame::create(uri, schema, index_columns, ctx);

    // Test NDRectangle functionality using SOMADataFrame
    {
        auto sdf = SOMADataFrame::open(uri, OpenMode::soma_read, ctx);
        auto raw_array = tiledb::Array(*ctx->tiledb_ctx(), uri, TILEDB_READ);
        auto columns = SOMAColumn::deserialize(*ctx->tiledb_ctx(), raw_array, {});

        // Find the string dimension column
        std::shared_ptr<SOMAColumn> string_dim_column = nullptr;
        for (const auto& column : columns) {
            if (column->tiledb_dimensions().has_value()) {
                auto dimensions = column->tiledb_dimensions().value();
                for (const auto& dimension : dimensions) {
                    if (dimension.name() == "str_dim") {
                        string_dim_column = column;
                        break;
                    }
                }
            }
        }

        REQUIRE(string_dim_column != nullptr);

        // Test current domain functionality
        auto current_domain_pair = string_dim_column->core_current_domain_slot<std::string>(*ctx, raw_array);
        REQUIRE(current_domain_pair.first == "");
        REQUIRE(current_domain_pair.second == "");

        // Test with NDRectangle if current domain exists
        if (sdf->has_current_domain()) {
            auto current_domain = sdf->get_current_domain_for_test();
            if (!current_domain.is_empty() && current_domain.type() == TILEDB_NDRECTANGLE) {
                NDRectangle ndrect = current_domain.ndrectangle();

                // Test NDRectangle version of the method
                auto ndrect_domain_pair = string_dim_column->core_current_domain_slot<std::string>(ndrect);
                REQUIRE(ndrect_domain_pair.first == "");
                REQUIRE(ndrect_domain_pair.second == "");
            }
        }

        sdf->close();
        raw_array.close();
    }
}
