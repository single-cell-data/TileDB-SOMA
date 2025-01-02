/**
 * @file   unit_soma_column.cc
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
 * This file manages unit tests for implementation of SOMAColumn class. This is
 * temparary and to be removed once SOMAColumn is fully integrated.
 */

#include <format>
#include <tiledb/tiledb>
#include <tiledbsoma/tiledbsoma>
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

    std::string i64_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        i64_datatype);
    std::string u32_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        u32_datatype);
    std::string attr_1_arrow_format = ArrowAdapter::tdb_to_arrow_type(
        str_datatype);

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
             {.name = "x",
              .tiledb_datatype = TILEDB_FLOAT64,
              .dim_max = 200,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "y",
              .tiledb_datatype = TILEDB_FLOAT64,
              .dim_max = 100,
              .string_lo = "N/A",
              .string_hi = "N/A"})});

    auto index_columns = helper::create_column_index_info(dim_infos);

    std::vector<std::shared_ptr<SOMAColumn>> columns;

    for (int64_t i = 0; i < index_columns.second->n_children; ++i) {
        columns.push_back(SOMADimension::create(
            ctx->tiledb_ctx(),
            index_columns.second->children[i],
            index_columns.first->children[i],
            "SOMAGeometryDataFrame",
            "",
            platform_config));

        REQUIRE(
            columns.back()->tiledb_dimensions().value()[0].type() ==
            dim_infos[i].tiledb_datatype);
    }

    REQUIRE(
        columns[1]->core_domain_slot<double_t>() ==
        std::make_pair<double_t, double_t>(0, helper::CORE_DOMAIN_MAX));
    REQUIRE(
        columns[1]->core_domain_slot<double_t>() ==
        std::make_pair<double_t, double_t>(0, helper::CORE_DOMAIN_MAX));
    REQUIRE(
        columns[2]->core_domain_slot<int64_t>() ==
        std::make_pair<int64_t, int64_t>(0, helper::CORE_DOMAIN_MAX));
    REQUIRE(
        columns[3]->core_domain_slot<std::string>() ==
        std::make_pair<std::string, std::string>("", ""));
}

TEST_CASE_METHOD(
    VariouslyIndexedDataFrameFixture,
    "SOMAColumn: query variant-indexed dataframe dim-str-u32 attr-sjid",
    "[SOMADataFrame]") {
    auto specify_domain = GENERATE(false, true);
    SECTION(std::format("- specify_domain={}", specify_domain)) {
        std::string suffix1 = specify_domain ? "true" : "false";
        set_up(
            std::make_shared<SOMAContext>(),
            "mem://unit-test-column-variant-indexed-dataframe-4-" + suffix1);

        std::string string_lo = "";
        std::string string_hi = "";
        std::vector<helper::DimInfo> dim_infos(
            {str_dim_info(string_lo, string_hi), u32_dim_info()});
        std::vector<helper::AttrInfo> attr_infos({i64_attr_info()});

        // Create
        create(dim_infos, attr_infos);

        // Check current domain
        auto sdf = open(OpenMode::read);

        // External column initialization
        auto raw_array = tiledb::Array(*ctx_->tiledb_ctx(), uri_, TILEDB_READ);
        std::vector<std::shared_ptr<SOMAColumn>> columns;

        for (auto dimension : sdf->tiledb_schema()->domain().dimensions()) {
            columns.push_back(
                std::make_shared<SOMADimension>(SOMADimension(dimension)));
        }

        CurrentDomain current_domain = sdf->get_current_domain_for_test();

        REQUIRE(!current_domain.is_empty());
        REQUIRE(current_domain.type() == TILEDB_NDRECTANGLE);
        NDRectangle ndrect = current_domain.ndrectangle();

        std::array<std::string, 2> str_range = ndrect.range<std::string>(
            dim_infos[0].name);
        std::pair<std::string, std::string>
            str_external = columns[0]->core_current_domain_slot<std::string>(
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

        std::array<uint32_t, 2> u32_range = ndrect.range<uint32_t>(
            dim_infos[1].name);
        std::pair<uint32_t, uint32_t>
            u32_external = columns[1]->core_current_domain_slot<uint32_t>(
                *ctx_, raw_array);
        REQUIRE(u32_range[0] == u32_external.first);
        REQUIRE(u32_range[1] == u32_external.second);

        // Check shape before write
        std::optional<int64_t> actual = sdf->maybe_soma_joinid_shape();
        REQUIRE(!actual.has_value());

        // Check domainish accessors before resize
        ArrowTable non_empty_domain = sdf->get_non_empty_domain();
        std::vector<std::string>
            ned_str = ArrowAdapter::get_table_string_column_by_name(
                non_empty_domain, "mystring");

        std::vector<std::string>
            ned_str_col = ArrowAdapter::get_array_string_column(
                columns[0]->arrow_domain_slot(
                    *ctx_, raw_array, Domainish::kind_non_empty_domain),
                columns[0]->arrow_schema_slot(*ctx_, raw_array));

        ArrowTable soma_domain = sdf->get_soma_domain();
        std::vector<std::string>
            dom_str = ArrowAdapter::get_table_string_column_by_name(
                soma_domain, "mystring");

        std::vector<std::string>
            dom_str_col = ArrowAdapter::get_array_string_column(
                columns[0]->arrow_domain_slot(
                    *ctx_, raw_array, Domainish::kind_core_current_domain),
                columns[0]->arrow_schema_slot(*ctx_, raw_array));

        ArrowTable soma_maxdomain = sdf->get_soma_maxdomain();
        std::vector<std::string>
            maxdom_str = ArrowAdapter::get_table_string_column_by_name(
                soma_maxdomain, "mystring");

        std::vector<std::string>
            maxdom_str_col = ArrowAdapter::get_array_string_column(
                columns[0]->arrow_domain_slot(
                    *ctx_, raw_array, Domainish::kind_core_domain),
                columns[0]->arrow_schema_slot(*ctx_, raw_array));

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

        sdf = open(OpenMode::write);
        write_sjid_u32_str_data_from(0);

        sdf->close();

        sdf = open(OpenMode::read);
        REQUIRE(sdf->nnz() == 2);

        sdf->close();

        auto external_query = std::make_unique<ManagedQuery>(
            open(OpenMode::read), ctx_->tiledb_ctx());

        columns[1]->select_columns(external_query);
        columns[1]->set_dim_point<uint32_t>(external_query, *ctx_, 1234);

        // Configure query and allocate result buffers
        auto ext_res = external_query->read_next().value();

        REQUIRE(ext_res->num_rows() == 1);

        external_query->reset();

        columns[0]->select_columns(external_query);
        columns[0]->set_dim_ranges<std::string>(
            external_query,
            *ctx_,
            std::vector(
                {std::make_pair<std::string, std::string>("apple", "b")}));

        // Configure query and allocate result buffers
        ext_res = external_query->read_next().value();

        REQUIRE(ext_res->num_rows() == 1);
    }
}
