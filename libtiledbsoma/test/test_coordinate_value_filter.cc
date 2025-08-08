/**
 * @file   test_coordinate_value_filter.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMAArray class
 */

#include "common.h"

/**
 * Test that checks SOMArdQueryCondition is grabbing the correct values.
 *
 * This is an ideal test to parametrize with rapidtest when we add it to SOMA.
 */
TEST_CASE("Test CoordinateValueFilters on SparseArray", "[CoordinateValueFilters][SOMASparseNDArray]") {
    // Create a TileDB sparse array with 3 integer dimensions.
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://test-query-condition-sparse";
    auto index_columns = helper::create_column_index_info(
        {helper::DimInfo(
             {.name = "soma_dim_0",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = 4,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "soma_dim_1",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = 3,
              .string_lo = "N/A",
              .string_hi = "N/A"})});

    SOMASparseNDArray::create(uri, "i", index_columns, ctx);
    std::vector<std::string> dim_names{"soma_dim_0", "soma_dim_1"};

    // Define input data.
    std::vector<int64_t> coords_dim_0{0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
    std::vector<int64_t> coords_dim_1{0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<int32_t> data(12);
    std::iota(data.begin(), data.end(), 1);

    // Write data to the array.
    {
        Array array{*ctx->tiledb_ctx(), uri, TILEDB_WRITE};
        Query query{*ctx->tiledb_ctx(), array};
        query.set_layout(TILEDB_GLOBAL_ORDER);
        query.set_data_buffer("soma_dim_0", coords_dim_0);
        query.set_data_buffer("soma_dim_1", coords_dim_1);
        query.set_data_buffer("soma_data", data);
        query.submit();
        query.finalize();
    }

    auto soma_array = SOMASparseNDArray::open(uri, OpenMode::soma_read, ctx);
    auto array = *soma_array->tiledb_array();
    auto check_query_condition =
        [&](const CoordinateValueFilters& value_filter, const Subarray& subarray, const std::string& log_note) {
            INFO(log_note);

            // Create a query for the entire array.
            std::vector<int64_t> expected_coords_0(12);
            std::vector<int64_t> expected_coords_1(12);
            std::vector<int32_t> expected_data(12);

            // Create a query for the entire array.
            std::vector<int64_t> actual_coords_0(12);
            std::vector<int64_t> actual_coords_1(12);
            std::vector<int32_t> actual_data(12);

            // Query with subarray.
            Query query1(*ctx->tiledb_ctx(), array);
            query1.set_layout(TILEDB_ROW_MAJOR);
            query1.set_data_buffer("soma_dim_0", expected_coords_0);
            query1.set_data_buffer("soma_dim_1", expected_coords_1);
            query1.set_data_buffer("soma_data", expected_data);
            query1.set_subarray(subarray);
            query1.submit();

            // Query with query condition.
            Query query2(*ctx->tiledb_ctx(), array);
            query2.set_layout(TILEDB_ROW_MAJOR);
            query2.set_data_buffer("soma_dim_0", actual_coords_0);
            query2.set_data_buffer("soma_dim_1", actual_coords_1);
            query2.set_data_buffer("soma_data", actual_data);
            Subarray value_filter_subarray(*ctx->tiledb_ctx(), array);
            value_filter_subarray.add_range<int64_t>(0, 0, 3).add_range<int64_t>(1, 0, 2);
            query2.set_subarray(value_filter_subarray);
            query2.set_condition(value_filter.combine().query_condition());
            query2.submit();

            // Check results.
            auto expected_result_num = static_cast<int64_t>(query1.result_buffer_elements()["soma_data"].second);
            auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_data"].second);
            CHECK(expected_result_num == actual_result_num);
            CHECK_THAT(actual_coords_0, Catch::Matchers::Equals(expected_coords_0));
            CHECK_THAT(actual_coords_1, Catch::Matchers::Equals(expected_coords_1));
            CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
        };

    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        Subarray subarray(*ctx->tiledb_ctx(), array);
        value_filter.add_column_selection<int64_t>(0, SliceSelection<int64_t>(0, 3))
            .add_column_selection<int64_t>(1, SliceSelection<int64_t>(0, 2));
        subarray.add_range<int64_t>(0, 0, 3).add_range<int64_t>(1, 0, 2);
        check_query_condition(value_filter, subarray, "Read all values by range.");
    }

    // Error: out-of-bounds range.
    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        CHECK_THROWS_AS(
            value_filter.add_column_selection<int64_t>(0, SliceSelection<int64_t>(5, 7)), std::out_of_range);
    }

    // Error: out-of-bounds points.
    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        std::vector<int64_t> points{5, 7, 11, 10};
        CHECK_THROWS_AS(
            value_filter.add_column_selection<int64_t>(0, PointSelection<int64_t>(points)), std::out_of_range);
    }

    {
        INFO("Error: slice with incorrect type");
        auto value_filter = soma_array->create_coordinate_value_filter();
        std::vector<int32_t> points{1, 4, 7};
        CHECK_THROWS_AS(value_filter.add_points<int32_t>(0, PointSelection<int32_t>(points)), std::invalid_argument);
    }
    {
        INFO("Error: points with incorrect type");
        auto value_filter = soma_array->create_coordinate_value_filter();
        CHECK_THROWS_AS(
            value_filter.add_column_selection<uint64_t>(0, SliceSelection<uint64_t>(0, 10)), std::invalid_argument);
    }

    // Region [1:2]x[:] by ranges.
    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        Subarray subarray(*ctx->tiledb_ctx(), array);
        value_filter.add_column_selection<int64_t>(0, SliceSelection<int64_t>(1, 2));
        subarray.add_range<int64_t>(0, 1, 2);
        check_query_condition(value_filter, subarray, "Select by range on dim 0.");
    }

    // Region [:]x[1:2] by ranges.
    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        Subarray subarray(*ctx->tiledb_ctx(), array);
        value_filter.add_column_selection<int64_t>(1, SliceSelection<int64_t>(1, 2));
        subarray.add_range<int64_t>(1, 1, 2);
        check_query_condition(value_filter, subarray, "Select by range on dim 1.");
    }

    // Region [0,1,3]x[:] by points (ordered).
    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        Subarray subarray(*ctx->tiledb_ctx(), array);
        std::vector<int64_t> points{0, 1, 3};
        value_filter.add_column_selection<int64_t>(0, PointSelection<int64_t>(points));
        subarray.add_range<int64_t>(0, 0, 1).add_range<int64_t>(0, 3, 3);
        check_query_condition(value_filter, subarray, "Select by points on dim 0 (ordered).");
    }
    // Region [0,1,3]x[:] by points (unordered).
    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        Subarray subarray(*ctx->tiledb_ctx(), array);
        std::vector<int64_t> points{3, 0, 1};
        value_filter.add_column_selection<int64_t>(0, PointSelection<int64_t>(points));
        subarray.add_range<int64_t>(0, 0, 1).add_range<int64_t>(0, 3, 3);
        check_query_condition(value_filter, subarray, "Select by points on dim 0 (unordered).");
    }

    // Region [0,1,3]x[:] by points (multiple conditions).
    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        Subarray subarray(*ctx->tiledb_ctx(), array);
        std::vector<int64_t> points1{3};
        std::vector<int64_t> points2{0};
        std::vector<int64_t> points3{1};
        value_filter.add_column_selection<int64_t>(0, PointSelection<int64_t>(points1))
            .add_column_selection<int64_t>(0, PointSelection<int64_t>(points2))
            .add_column_selection<int64_t>(0, PointSelection<int64_t>(points3));
        subarray.add_range<int64_t>(0, 0, 1).add_range<int64_t>(0, 3, 3);
        check_query_condition(value_filter, subarray, "Select by points on dim 0 (multiple conditions).");
    }

    // Region [:]x[0,2] by points.
    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        Subarray subarray(*ctx->tiledb_ctx(), array);
        std::vector<int64_t> points{0, 2};
        value_filter.add_column_selection<int64_t>(1, PointSelection<int64_t>(points));
        subarray.add_range<int64_t>(1, 0, 0).add_range<int64_t>(1, 2, 2);
        check_query_condition(value_filter, subarray, "Select by points on dim 1.");
    }
}

TEST_CASE(
    "Test CoordinateValueFilters on SOMADataFrame with string index column",
    "[CoordinateValueFilters][SOMADataFrame][string-index]") {
    // Create a TileDB sparse array with 1 integer dimensions.
    auto ctx = std::make_shared<SOMAContext>();
    auto tiledb_ctx = ctx->tiledb_ctx();
    std::string uri = "mem://test-query-condition-string";
    auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(
        {helper::DimInfo(
            {.name = "label", .tiledb_datatype = TILEDB_STRING_ASCII, .dim_max = 0, .string_lo = "", .string_hi = ""})},
        {helper::AttrInfo({.name = "soma_joinid", .tiledb_datatype = TILEDB_INT64})});

    SOMADataFrame::create(uri, schema, index_columns, ctx);

    // Define input data.
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

    // Write data to the array.
    {
        Array array{*ctx->tiledb_ctx(), uri, TILEDB_WRITE};
        Query query{*ctx->tiledb_ctx(), array};
        query.set_data_buffer("label", coords_data);
        query.set_offsets_buffer("label", coords_offsets);
        query.set_data_buffer("soma_joinid", index);
        query.submit();
        query.finalize();
    }

    auto soma_array = SOMADataFrame::open(uri, OpenMode::soma_read, ctx);
    auto array = *soma_array->tiledb_array();
    auto check_query_condition =
        [&](const CoordinateValueFilters& value_filter, const Subarray& subarray, const std::string& log_note) {
            INFO(log_note);

            // Check valid query condition.
            REQUIRE(value_filter.is_initialized());

            // Create a query for the entire array.
            std::vector<char> expected_coord_data(35);
            std::vector<uint64_t> expected_coord_offsets(7);
            std::vector<int64_t> expected_index(6);

            std::vector<char> actual_coord_data(35);
            std::vector<uint64_t> actual_coord_offsets(7);
            std::vector<int64_t> actual_index(6);

            // Query with subarray.
            Query query1(*ctx->tiledb_ctx(), array);
            query1.set_layout(TILEDB_ROW_MAJOR);
            query1.set_data_buffer("label", expected_coord_data);
            query1.set_offsets_buffer("label", expected_coord_offsets);
            query1.set_data_buffer("soma_joinid", expected_index);
            query1.set_subarray(subarray);
            query1.submit();
            REQUIRE(query1.query_status() == tiledb::Query::Status::COMPLETE);

            // Query with query condition.
            Query query2(*ctx->tiledb_ctx(), array);
            query2.set_layout(TILEDB_ROW_MAJOR);
            query2.set_data_buffer("label", actual_coord_data);
            query2.set_offsets_buffer("label", actual_coord_offsets);
            query2.set_data_buffer("soma_joinid", actual_index);
            Subarray value_filter_subarray(*ctx->tiledb_ctx(), array);
            value_filter_subarray.add_range(0, std::string("a"), std::string("z"));
            query2.set_subarray(value_filter_subarray);
            query2.set_condition(value_filter.combine().query_condition());
            query2.submit();
            REQUIRE(query2.query_status() == tiledb::Query::Status::COMPLETE);

            // Check results.
            auto expected_result_num = static_cast<int64_t>(query1.result_buffer_elements()["soma_joinid"].second);
            auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_joinid"].second);
            CHECK(expected_result_num == actual_result_num);
            CHECK_THAT(actual_coord_data, Catch::Matchers::Equals(expected_coord_data));
            CHECK_THAT(actual_coord_offsets, Catch::Matchers::Equals(expected_coord_offsets));
            CHECK_THAT(actual_index, Catch::Matchers::Equals(expected_index));
        };

    auto check_empty_query_condition = [&](const CoordinateValueFilters& value_filter, const std::string& log_note) {
        INFO(log_note);

        // Check valid query condition.
        REQUIRE(value_filter.is_initialized());

        // Create a query for the entire array.
        std::vector<char> expected_coord_data(35);
        std::vector<uint64_t> expected_coord_offsets(7);
        std::vector<int64_t> expected_index(6);

        std::vector<char> actual_coord_data(35);
        std::vector<uint64_t> actual_coord_offsets(7);
        std::vector<int64_t> actual_index(6);

        // Query with query condition.
        Query query2(*ctx->tiledb_ctx(), array);
        query2.set_layout(TILEDB_ROW_MAJOR);
        query2.set_data_buffer("label", actual_coord_data);
        query2.set_offsets_buffer("label", actual_coord_offsets);
        query2.set_data_buffer("soma_joinid", actual_index);
        Subarray value_filter_subarray(*ctx->tiledb_ctx(), array);
        value_filter_subarray.add_range(0, std::string("a"), std::string("z"));
        query2.set_subarray(value_filter_subarray);
        query2.set_condition(value_filter.combine().query_condition());
        query2.submit();
        REQUIRE(query2.query_status() == tiledb::Query::Status::COMPLETE);

        // Check results.
        auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_joinid"].second);
        CHECK(actual_result_num == 0);
        CHECK_THAT(actual_coord_data, Catch::Matchers::Equals(expected_coord_data));
        CHECK_THAT(actual_coord_offsets, Catch::Matchers::Equals(expected_coord_offsets));
        CHECK_THAT(actual_index, Catch::Matchers::Equals(expected_index));
    };
    std::vector<std::string> dim_names{"label"};

    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        value_filter.add_column_selection<std::string>(0, SliceSelection<std::string>("a", "z"));
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("apple"), std::string("fig"));
        check_query_condition(value_filter, subarray, "Read all values by range.");
    }

    {
        INFO("Invalid points: no values selected - should throw an std::invalid_argument error");
        auto value_filter = soma_array->create_coordinate_value_filter();
        std::vector<std::string> points{};
        CHECK_THROWS_AS(
            value_filter.add_column_selection<std::string>(0, PointSelection<std::string>(points)),
            std::invalid_argument);
    }
    {
        INFO("Error: slice with incorrect type");
        auto value_filter = soma_array->create_coordinate_value_filter();
        std::vector<int32_t> points{1, 4, 7};
        CHECK_THROWS_AS(
            value_filter.add_column_selection<int32_t>(0, PointSelection<int32_t>(points)), std::invalid_argument);
    }
    {
        INFO("Error: points with incorrect type");
        auto value_filter = soma_array->create_coordinate_value_filter();
        CHECK_THROWS_AS(
            value_filter.add_column_selection<int64_t>(0, SliceSelection<int64_t>(-1.5, std::nullopt)),
            std::invalid_argument);
    }

    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        value_filter.add_column_selection<std::string>(0, SliceSelection<std::string>("ca", "fa"));
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("ca"), std::string("fa"));
        check_query_condition(value_filter, subarray, "Create by range.");
    }

    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        std::vector<std::string> points{"fig", "durian", "banana"};
        value_filter.add_column_selection<std::string>(0, PointSelection<std::string>(points));
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("banana"), std::string("banana"))
            .add_range(0, std::string("durian"), std::string("durian"))
            .add_range(0, std::string("fig"), std::string("fig"));
        check_query_condition(value_filter, subarray, "Select by points.");
    }

    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        value_filter.add_column_selection<std::string>(0, SliceSelection<std::string>("g", "z"));
        check_empty_query_condition(value_filter, "Create by range - no values selected.");
    }

    {
        auto value_filter = soma_array->create_coordinate_value_filter();
        std::vector<std::string> points{"kiwi", "pear", "carrot"};
        value_filter.add_column_selection<std::string>(0, PointSelection<std::string>(points));
        check_empty_query_condition(value_filter, "Create by points - no values selected.");
    }
}
