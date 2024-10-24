/**
 * @file   unit_soma_geometry_dataframe.cc
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
 * This file manages unit tests for the SOMAGeometryDataFrame class
 */

#include <vector>
#include "../src/geometry/geometry.h"
#include "../src/geometry/operators/io/write.h"
#include "../src/utils/common.h"
#include "common.h"

const int64_t SOMA_JOINID_DIM_MAX = 99;

TEST_CASE("SOMAGeometryDataFrame: basic", "[SOMAGeometryDataFrame]") {
    auto use_current_domain = GENERATE(false, true);
    // TODO this could be formatted with fmt::format which is part of internal
    // header spd/log/fmt/fmt.h and should not be used. In C++20, this can be
    // replaced with std::format.
    std::ostringstream section;
    section << "- use_current_domain=" << use_current_domain;
    SECTION(section.str()) {
        auto ctx = std::make_shared<SOMAContext>();
        std::string uri{"mem://unit-test-geometry-basic"};
        PlatformConfig platform_config{};

        std::vector<helper::DimInfo> dim_infos(
            {helper::DimInfo(
                 {.name = "soma_joinid",
                  .tiledb_datatype = TILEDB_INT64,
                  .dim_max = SOMA_JOINID_DIM_MAX,
                  .string_lo = "N/A",
                  .string_hi = "N/A",
                  .use_current_domain = use_current_domain}),
             helper::DimInfo(
                 {.name = "soma_geometry",
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

        std::vector<helper::AttrInfo> attr_infos({helper::AttrInfo(
            {.name = "quality", .tiledb_datatype = TILEDB_FLOAT64})});

        // Check the geometry dataframe doesn't exist yet.
        REQUIRE(!SOMAGeometryDataFrame::exists(uri, ctx));

        // Create the geometry dataframe.
        auto [schema, index_columns] =
            helper::create_arrow_schema_and_index_columns(
                dim_infos, attr_infos);
        auto spatial_columns = helper::create_column_index_info(
            spatial_dim_infos);

        SOMAGeometryDataFrame::create(
            uri,
            std::move(schema),
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ArrowTable(
                std::move(spatial_columns.first),
                std::move(spatial_columns.second)),
            ctx,
            platform_config,
            std::nullopt);

        // Check the geometry dataframe exists and it cannot be read as a
        // different object.
        REQUIRE(SOMAGeometryDataFrame::exists(uri, ctx));
        REQUIRE(!SOMASparseNDArray::exists(uri, ctx));
        REQUIRE(!SOMADenseNDArray::exists(uri, ctx));
        REQUIRE(!SOMADataFrame::exists(uri, ctx));

        auto soma_geometry = SOMAGeometryDataFrame::open(
            uri,
            OpenMode::read,
            ctx,
            {},  // column_names,
            ResultOrder::automatic,
            std::nullopt);
        REQUIRE(soma_geometry->uri() == uri);
        REQUIRE(soma_geometry->ctx() == ctx);
        REQUIRE(soma_geometry->type() == "SOMAGeometryDataFrame");
        std::vector<std::string> expected_index_column_names = {
            dim_infos[0].name,
            SOMA_GEOMETRY_DIMENSION_PREFIX + spatial_dim_infos[0].name +
                "__min",
            SOMA_GEOMETRY_DIMENSION_PREFIX + spatial_dim_infos[1].name +
                "__min",
            SOMA_GEOMETRY_DIMENSION_PREFIX + spatial_dim_infos[0].name +
                "__max",
            SOMA_GEOMETRY_DIMENSION_PREFIX + spatial_dim_infos[1].name +
                "__max"};

        std::vector<std::string> expected_spatial_column_names = {
            spatial_dim_infos[0].name, spatial_dim_infos[1].name};
        REQUIRE(
            soma_geometry->index_column_names() == expected_index_column_names);
        REQUIRE(
            soma_geometry->spatial_column_names() ==
            expected_spatial_column_names);
        REQUIRE(soma_geometry->nnz() == 0);
        soma_geometry->close();

        auto soma_object = SOMAObject::open(uri, OpenMode::read, ctx);
        REQUIRE(soma_object->uri() == uri);
        REQUIRE(soma_object->type() == "SOMAGeometryDataFrame");
        soma_object->close();
    }
}