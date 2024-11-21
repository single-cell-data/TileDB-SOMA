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
 * This file manages unit tests for implementation of SOMAColumn class
 */

#include <tiledbsoma/tiledbsoma>
#include "../src/soma/soma_attribute.h"
#include "../src/soma/soma_column.h"
#include "../src/soma/soma_dimension.h"
#include "../src/soma/soma_geometry_column.h"
#include "common.h"

TEST_CASE("SOMAColumn: SOMADimension") {
    auto use_current_domain = GENERATE(false, true);
    auto ctx = std::make_shared<SOMAContext>();
    PlatformConfig platform_config{};

    SECTION(std::format("- use_current_domain={}", use_current_domain)) {
        std::vector<helper::DimInfo> dim_infos(
            {helper::DimInfo(
                 {.name = "dimension",
                  .tiledb_datatype = TILEDB_UINT32,
                  .dim_max = 100,
                  .string_lo = "N/A",
                  .string_hi = "N/A",
                  .use_current_domain = use_current_domain}),
             helper::DimInfo(
                 {.name = "dimension",
                  .tiledb_datatype = TILEDB_FLOAT64,
                  .dim_max = 100,
                  .string_lo = "N/A",
                  .string_hi = "N/A",
                  .use_current_domain = use_current_domain}),
             helper::DimInfo(
                 {.name = "dimension",
                  .tiledb_datatype = TILEDB_INT64,
                  .dim_max = 100,
                  .string_lo = "N/A",
                  .string_hi = "N/A",
                  .use_current_domain = use_current_domain}),
             helper::DimInfo(
                 {.name = "dimension",
                  .tiledb_datatype = TILEDB_STRING_ASCII,
                  .dim_max = 100,
                  .string_lo = "N/A",
                  .string_hi = "N/A",
                  .use_current_domain = use_current_domain})});

        std::vector<helper::DimInfo> geom_dim_infos({helper::DimInfo(
            {.name = "dimension",
             .tiledb_datatype = TILEDB_GEOM_WKB,
             .dim_max = 100,
             .string_lo = "N/A",
             .string_hi = "N/A",
             .use_current_domain = use_current_domain})});

        std::vector<helper::DimInfo> spatial_dim_infos(
            {helper::DimInfo(
                 {.name = "x",
                  .tiledb_datatype = TILEDB_FLOAT64,
                  .dim_max = 200,
                  .string_lo = "N/A",
                  .string_hi = "N/A",
                  .use_current_domain = use_current_domain}),
             helper::DimInfo(
                 {.name = "y",
                  .tiledb_datatype = TILEDB_FLOAT64,
                  .dim_max = 100,
                  .string_lo = "N/A",
                  .string_hi = "N/A",
                  .use_current_domain = use_current_domain})});

        auto index_columns = helper::create_column_index_info(dim_infos);

        std::vector<std::shared_ptr<SOMAColumn>> columns;
        bool has_current_domain = true;

        for (int64_t i = 0; i < index_columns.second->n_children; ++i) {
            columns.push_back(SOMADimension::create(
                ctx->tiledb_ctx(),
                index_columns.second->children[i],
                index_columns.first->children[i],
                "SOMAGeometryDataFrame",
                "",
                platform_config,
                has_current_domain));

            REQUIRE(has_current_domain == use_current_domain);
            REQUIRE(
                columns.back()->tiledb_dimensions().value()[0].type() ==
                dim_infos[i].tiledb_datatype);
        }

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
    "[SOMAColumn]") {
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

        for (size_t i = 0; i < sdf->tiledb_schema()->attribute_num(); ++i) {
            columns.push_back(std::make_shared<SOMAAttribute>(
                SOMAAttribute(sdf->tiledb_schema()->attribute(i))));
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
