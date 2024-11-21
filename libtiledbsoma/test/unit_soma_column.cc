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

#include <tiledbsoma/tiledbsoma>
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
            columns[1]->core_domain_slot<double_t>() ==
            std::make_pair<double_t, double_t>(
                0, use_current_domain ? helper::CORE_DOMAIN_MAX : 100));
        REQUIRE(
            columns[1]->core_domain_slot<double_t>() ==
            std::make_pair<double_t, double_t>(
                0, use_current_domain ? helper::CORE_DOMAIN_MAX : 100));
        REQUIRE(
            columns[2]->core_domain_slot<int64_t>() ==
            std::make_pair<int64_t, int64_t>(
                0, use_current_domain ? helper::CORE_DOMAIN_MAX : 100));
        REQUIRE(
            columns[3]->core_domain_slot<std::string>() ==
            std::make_pair<std::string, std::string>("", ""));
    }
}

TEST_CASE("SOMAColumn: SOMAGeometryDimension") {
    auto use_current_domain = GENERATE(false, true);
    auto ctx = std::make_shared<SOMAContext>();
    PlatformConfig platform_config{};

    SECTION(std::format("- use_current_domain={}", use_current_domain)) {
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

        auto geom_columns = helper::create_column_index_info(geom_dim_infos);
        auto spatial_columns = helper::create_column_index_info(
            spatial_dim_infos);

        bool has_current_domain = true;
        auto geometry_column = SOMAGeometryColumn::create(
            ctx->tiledb_ctx(),
            geom_columns.second->children[0],
            geom_columns.first->children[0],
            spatial_columns.second.get(),
            spatial_columns.first.get(),
            "SOMAGeometryDataFrame",
            "WKB",
            platform_config,
            has_current_domain);

        REQUIRE(
            geometry_column->tiledb_dimensions().value().size() ==
            spatial_dim_infos.size() * 2);
        REQUIRE(
            geometry_column->tiledb_attributes().value()[0].type() ==
            TILEDB_GEOM_WKB);

        auto domain = geometry_column
                          ->core_domain_slot<std::vector<double_t>>();
        CHECK(domain.first == std::vector<double_t>({0, 0}));
        CHECK(
            domain.second ==
            std::vector<double_t>(
                {use_current_domain ? helper::CORE_DOMAIN_MAX : 200.0,
                 use_current_domain ? helper::CORE_DOMAIN_MAX : 100.0}));
    }
}