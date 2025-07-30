/**
 * @file   unit_soma_query_condition.cc
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
#include "tiledb_adapter/soma_query_condition.h"

/**
 * Test that checks SOMAValueFilter is grabbing the correct values.
 *
 * This is an ideal test to parametrize with rapidtest when we add it to SOMA.
 */
TEST_CASE("Test SOMAValueFilter on SparseArray", "[SOMAValueFilter][SOMASparseNDArray]") {
    // Create a TileDB sparse array with 1 integer dimensions.
    auto ctx = std::make_shared<SOMAContext>();
    auto tiledb_ctx = ctx->tiledb_ctx();
    std::string uri = "mem://test-query-condition-sparse";
    auto index_columns = helper::create_column_index_info({helper::DimInfo(
        {.name = "soma_dim_0",
         .tiledb_datatype = TILEDB_INT64,
         .dim_max = 16,
         .string_lo = "N/A",
         .string_hi = "N/A"})});

    SOMASparseNDArray::create(uri, "i", index_columns, ctx);

    // Define input data.
    std::vector<int64_t> coords(16);
    std::vector<int32_t> data(16);
    std::iota(coords.begin(), coords.end(), 0);
    std::iota(data.begin(), data.end(), 1);

    // Write data to the array.
    {
        Array array{*ctx->tiledb_ctx(), uri, TILEDB_WRITE};
        Query query{*ctx->tiledb_ctx(), array};
        query.set_layout(TILEDB_GLOBAL_ORDER);
        query.set_data_buffer("soma_dim_0", coords);
        query.set_data_buffer("soma_data", data);
        query.submit();
        query.finalize();
    }

    Array array{*ctx->tiledb_ctx(), uri, TILEDB_READ};
    auto check_query_condition = [&](const SOMAValueFilter& qc, const Subarray& subarray, const std::string& log_note) {
        INFO(log_note);

        // Create a query for the entire array.
        std::vector<int64_t> expected_coords_0(16);
        std::vector<int32_t> expected_data(16);

        // Create a query for the entire array.
        std::vector<int64_t> actual_coords_0(16);
        std::vector<int32_t> actual_data(16);

        // Query with subarray.
        Query query1(*ctx->tiledb_ctx(), array);
        query1.set_layout(TILEDB_ROW_MAJOR);
        query1.set_data_buffer("soma_dim_0", expected_coords_0);
        query1.set_data_buffer("soma_data", expected_data);
        query1.set_subarray(subarray);
        query1.submit();

        // Query with query condition.
        Query query2(*ctx->tiledb_ctx(), array);
        query2.set_layout(TILEDB_ROW_MAJOR);
        query2.set_data_buffer("soma_dim_0", actual_coords_0);
        query2.set_data_buffer("soma_data", actual_data);
        Subarray qc_subarray(*ctx->tiledb_ctx(), array);
        qc_subarray.add_range<int64_t>(0, 0, 15);
        query2.set_subarray(qc_subarray);
        query2.set_condition(qc.query_condition());
        query2.submit();

        // Check results.
        auto expected_result_num = static_cast<int64_t>(query1.result_buffer_elements()["soma_data"].second);
        auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_data"].second);
        CHECK(expected_result_num == actual_result_num);
        CHECK_THAT(actual_coords_0, Catch::Matchers::Equals(expected_coords_0));
        CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
    };

    auto check_empty_query_condition = [&](const SOMAValueFilter& qc, const std::string& log_note) {
        INFO(log_note);

        // Create a query for the entire array.
        std::vector<int64_t> expected_coords_0(16);
        std::vector<int32_t> expected_data(16);

        // Create a query for the entire array.
        std::vector<int64_t> actual_coords_0(16);
        std::vector<int32_t> actual_data(16);

        // Query with query condition.
        Query query2(*ctx->tiledb_ctx(), array);
        query2.set_layout(TILEDB_ROW_MAJOR);
        query2.set_data_buffer("soma_dim_0", actual_coords_0);
        query2.set_data_buffer("soma_data", actual_data);
        Subarray qc_subarray(*ctx->tiledb_ctx(), array);
        qc_subarray.add_range<int64_t>(0, 0, 15);
        query2.set_subarray(qc_subarray);
        query2.set_condition(qc.query_condition());
        query2.submit();

        // Check results.
        auto expected_result_num = 0;
        auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_data"].second);
        CHECK(expected_result_num == actual_result_num);
        CHECK_THAT(actual_coords_0, Catch::Matchers::Equals(expected_coords_0));
        CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
    };

    // Check default constructor.
    {
        SOMAValueFilter qc{};
        CHECK(!qc.is_initialized());
    }

    // Full region by range.
    {
        auto qc = SOMAValueFilter::create_from_slice<int64_t>(*tiledb_ctx, "soma_dim_0", 0, 15);
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 0, 15);
        check_query_condition(qc, subarray, "Read all values by range.");
    }

    // Error: invalid range.
    {
        CHECK_THROWS_AS(
            SOMAValueFilter::create_from_slice<int64_t>(*tiledb_ctx, "soma_dim_0", 3, 2), std::invalid_argument);
    }

    // Empty region: out-of-bounds range.
    {
        auto qc = SOMAValueFilter::create_from_slice<int64_t>(*tiledb_ctx, "soma_dim_0", 25, 27);
        CHECK(qc.is_initialized());
        check_empty_query_condition(qc, "Range out of bounds: expect no values.");
    }

    // Empty region: invalid - no points.
    {
        std::vector<int64_t> points{};
        CHECK_THROWS_AS(
            SOMAValueFilter::create_from_points<int64_t>(*tiledb_ctx, "soma_dim_0", points), std::invalid_argument);
    }

    // Empty region: out-of-bounds points.
    {
        std::vector<int64_t> points{25, 27, 19, 30};
        auto qc = SOMAValueFilter::create_from_points<int64_t>(*tiledb_ctx, "soma_dim_0", points);
        CHECK(qc.is_initialized());
        check_empty_query_condition(qc, "Range out of bounds: expect no values.");
    }

    // Region [1:2] by ranges.
    {
        auto qc = SOMAValueFilter::create_from_slice<int64_t>(*tiledb_ctx, "soma_dim_0", 1, 2);
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 1, 2);
        check_query_condition(qc, subarray, "Select by range on dim 0.");
    }

    // Region [0,11,13] by points.
    {
        std::vector<int64_t> points{0, 1, 13};
        auto qc = SOMAValueFilter::create_from_points<int64_t>(*tiledb_ctx, "soma_dim_0", points);
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 0, 1).add_range<int64_t>("soma_dim_0", 13, 13);
        check_query_condition(qc, subarray, "Select by points on dim 0.");
    }
}

TEST_CASE(
    "Test SOMAValueFilter on SOMADataFrame with string index column",
    "[SOMAValueFilter][SOMADataFrame][string-index]") {
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

    Array array{*ctx->tiledb_ctx(), uri, TILEDB_READ};
    auto check_query_condition = [&](const SOMAValueFilter& qc, const Subarray& subarray, const std::string& log_note) {
        INFO(log_note);

        // Check valid query condition.
        REQUIRE(qc.is_initialized());

        // Create a query for the entire array.
        std::vector<char> expected_coord_data(data_size);
        std::vector<uint64_t> expected_coord_offsets(7);
        std::vector<int64_t> expected_index(6);

        std::vector<char> actual_coord_data(data_size);
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
        Subarray qc_subarray(*ctx->tiledb_ctx(), array);
        qc_subarray.add_range(0, std::string("a"), std::string("z"));
        query2.set_subarray(qc_subarray);
        query2.set_condition(qc.query_condition());
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

    auto check_empty_query_condition = [&](const SOMAValueFilter& qc, const std::string& log_note) {
        INFO(log_note);

        // Check valid query condition.
        REQUIRE(qc.is_initialized());

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
        Subarray qc_subarray(*ctx->tiledb_ctx(), array);
        qc_subarray.add_range(0, std::string("a"), std::string("z"));
        query2.set_subarray(qc_subarray);
        query2.set_condition(qc.query_condition());
        query2.submit();
        REQUIRE(query2.query_status() == tiledb::Query::Status::COMPLETE);

        // Check results.
        auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_joinid"].second);
        CHECK(actual_result_num == 0);
        CHECK_THAT(actual_coord_data, Catch::Matchers::Equals(expected_coord_data));
        CHECK_THAT(actual_coord_offsets, Catch::Matchers::Equals(expected_coord_offsets));
        CHECK_THAT(actual_index, Catch::Matchers::Equals(expected_index));
    };

    {
        auto qc = SOMAValueFilter::create_from_slice(*tiledb_ctx, "label", std::string("a"), std::string("z"));
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("apple"), std::string("fig"));
        check_query_condition(qc, subarray, "Read all values by range.");
    }

    {
        INFO("Invalid range: should throw an std::invalid_argument error");
        CHECK_THROWS_AS(
            SOMAValueFilter::create_from_slice(*tiledb_ctx, "label", std::string("fig"), std::string("apple")),
            std::invalid_argument);
    }

    {
        INFO("Invalid points: no values selected - should throw an std::invalid_argument error");
        std::vector<std::string> points{};
        CHECK_THROWS_AS(SOMAValueFilter::create_from_points(*tiledb_ctx, "label", points), std::invalid_argument);
    }

    {
        auto qc = SOMAValueFilter::create_from_slice(*tiledb_ctx, "label", std::string("ca"), std::string("fa"));
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("ca"), std::string("fa"));
        check_query_condition(qc, subarray, "Create by range.");
    }

    {
        std::vector<std::string> points{"fig", "durian", "banana"};
        auto qc = SOMAValueFilter::create_from_points(*tiledb_ctx, "label", points);
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("banana"), std::string("banana"))
            .add_range(0, std::string("durian"), std::string("durian"))
            .add_range(0, std::string("fig"), std::string("fig"));
        check_query_condition(qc, subarray, "Select by points.");
    }

    {
        auto qc = SOMAValueFilter::create_from_slice(*tiledb_ctx, "label", std::string("g"), std::string("z"));
        check_empty_query_condition(qc, "Create by range - no values selected.");
    }

    {
        std::vector<std::string> points{"kiwi", "pear", "carrot"};
        auto qc = SOMAValueFilter::create_from_points(*tiledb_ctx, "label", points);
        check_empty_query_condition(qc, "Create by points - no values selected.");
    }
}

/**
 * Test that checks SOMArdQueryCondition is grabbing the correct values.
 *
 * This is an ideal test to parametrize with rapidtest when we add it to SOMA.
 */
TEST_CASE("Test CoordinateValueFilter on SparseArray", "[CoordinateValueFilter][SOMASparseNDArray]") {
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

    Array array{*ctx->tiledb_ctx(), uri, TILEDB_READ};
    auto check_query_condition =
        [&](const CoordinateValueFilter& qc, const Subarray& subarray, const std::string& log_note) {
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
            Subarray qc_subarray(*ctx->tiledb_ctx(), array);
            qc_subarray.add_range<int64_t>(0, 0, 3).add_range<int64_t>(1, 0, 2);
            query2.set_subarray(qc_subarray);
            query2.set_condition(qc.get_value_filter().query_condition());
            query2.submit();

            // Check results.
            auto expected_result_num = static_cast<int64_t>(query1.result_buffer_elements()["soma_data"].second);
            auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_data"].second);
            CHECK(expected_result_num == actual_result_num);
            CHECK_THAT(actual_coords_0, Catch::Matchers::Equals(expected_coords_0));
            CHECK_THAT(actual_coords_1, Catch::Matchers::Equals(expected_coords_1));
            CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
        };

    auto check_empty_query_condition = [&](const CoordinateValueFilter& qc, const std::string& log_note) {
        INFO(log_note);

        // Create a query for the entire array.
        std::vector<int64_t> expected_coords_0(12);
        std::vector<int64_t> expected_coords_1(12);
        std::vector<int32_t> expected_data(12);

        // Create a query for the entire array.
        std::vector<int64_t> actual_coords_0(12);
        std::vector<int64_t> actual_coords_1(12);
        std::vector<int32_t> actual_data(12);

        // Query with query condition.
        Query query2(*ctx->tiledb_ctx(), array);
        query2.set_layout(TILEDB_ROW_MAJOR);
        query2.set_data_buffer("soma_dim_0", actual_coords_0);
        query2.set_data_buffer("soma_dim_1", actual_coords_1);
        query2.set_data_buffer("soma_data", actual_data);
        Subarray qc_subarray(*ctx->tiledb_ctx(), array);
        qc_subarray.add_range<int64_t>(0, 0, 3).add_range<int64_t>(1, 0, 2);
        query2.set_subarray(qc_subarray);
        query2.set_condition(qc.get_value_filter().query_condition());
        query2.submit();

        // Check results.
        auto expected_result_num = 0;
        auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_data"].second);
        CHECK(expected_result_num == actual_result_num);
        CHECK_THAT(actual_coords_0, Catch::Matchers::Equals(expected_coords_0));
        CHECK_THAT(actual_coords_1, Catch::Matchers::Equals(expected_coords_1));
        CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
    };

    // Full region by range.
    std::pair<int64_t, int64_t> domain0{0, 3};
    std::pair<int64_t, int64_t> domain1{0, 2};
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_column_selection<int64_t>(0, SOMASliceSelection<int64_t>(0, 3), domain0)
            .add_column_selection<int64_t>(1, SOMASliceSelection<int64_t>(0, 2), domain1);
        subarray.add_range<int64_t>(0, 0, 3).add_range<int64_t>(1, 0, 2);
        check_query_condition(qc, subarray, "Read all values by range.");
    }

    // Empty region: out-of-bounds range - incorrect domain provided
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        qc.add_column_selection<int64_t>(0, SOMASliceSelection<int64_t>(5, 7), std::pair<int64_t, int64_t>(0, 100));
        check_empty_query_condition(qc, "Range out of bounds: expect no values.");
    }

    // Error: out-of-bounds range.
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        CHECK_THROWS_AS(
            qc.add_column_selection<int64_t>(0, SOMASliceSelection<int64_t>(5, 7), domain0), std::out_of_range);
    }

    // Empty region: out-of-bounds points - incorrect domain provided
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        std::vector<int64_t> points{5, 7, 11, 10};
        qc.add_column_selection<int64_t>(0, SOMAPointSelection<int64_t>(points), std::pair<int64_t, int64_t>(0, 100));
        check_empty_query_condition(qc, "Points out of bounds: expect no values.");
    }

    // Error: out-of-bounds points.
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        std::vector<int64_t> points{5, 7, 11, 10};
        CHECK_THROWS_AS(
            qc.add_column_selection<int64_t>(0, SOMAPointSelection<int64_t>(points), domain0), std::out_of_range);
    }

    // Region [1:2]x[:] by ranges.
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_column_selection<int64_t>(0, SOMASliceSelection<int64_t>(1, 2), {0, 3});
        subarray.add_range<int64_t>(0, 1, 2);
        check_query_condition(qc, subarray, "Select by range on dim 0.");
    }

    // Region [:]x[1:2] by ranges.
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_column_selection<int64_t>(1, SOMASliceSelection<int64_t>(1, 2), {0, 2});
        subarray.add_range<int64_t>(1, 1, 2);
        check_query_condition(qc, subarray, "Select by range on dim 1.");
    }

    // Region [0,1,3]x[:] by points (ordered).
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        std::vector<int64_t> points{0, 1, 3};
        qc.add_column_selection<int64_t>(0, SOMAPointSelection<int64_t>(points), {0, 3});
        subarray.add_range<int64_t>(0, 0, 1).add_range<int64_t>(0, 3, 3);
        check_query_condition(qc, subarray, "Select by points on dim 0 (ordered).");
    }
    // Region [0,1,3]x[:] by points (unordered).
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        std::vector<int64_t> points{3, 0, 1};
        qc.add_column_selection<int64_t>(0, SOMAPointSelection<int64_t>(points), {0, 3});
        subarray.add_range<int64_t>(0, 0, 1).add_range<int64_t>(0, 3, 3);
        check_query_condition(qc, subarray, "Select by points on dim 0 (unordered).");
    }

    // Region [0,1,3]x[:] by points (multiple conditions).
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        std::vector<int64_t> points1{3};
        std::vector<int64_t> points2{0};
        std::vector<int64_t> points3{1};
        qc.add_column_selection<int64_t>(0, SOMAPointSelection<int64_t>(points1), {0, 3})
            .add_column_selection<int64_t>(0, SOMAPointSelection<int64_t>(points2), {0, 3})
            .add_column_selection<int64_t>(0, SOMAPointSelection<int64_t>(points3), {0, 3});
        subarray.add_range<int64_t>(0, 0, 1).add_range<int64_t>(0, 3, 3);
        check_query_condition(qc, subarray, "Select by points on dim 0 (multiple conditions).");
    }

    // Region [:]x[0,2] by points.
    {
        CoordinateValueFilter qc(*ctx, dim_names);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        std::vector<int64_t> points{0, 2};
        qc.add_column_selection<int64_t>(1, SOMAPointSelection<int64_t>(points), {0, 3});
        subarray.add_range<int64_t>(1, 0, 0).add_range<int64_t>(1, 2, 2);
        check_query_condition(qc, subarray, "Select by points on dim 1.");
    }
}

TEST_CASE(
    "Test CoordinateValueFilter on SOMADataFrame with string index column",
    "[SOMAValueFilter][SOMADataFrame][string-index]") {
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

    Array array{*ctx->tiledb_ctx(), uri, TILEDB_READ};
    auto check_query_condition =
        [&](const CoordinateValueFilter& qc, const Subarray& subarray, const std::string& log_note) {
            INFO(log_note);

            // Check valid query condition.
            REQUIRE(qc.is_initialized());

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
            Subarray qc_subarray(*ctx->tiledb_ctx(), array);
            qc_subarray.add_range(0, std::string("a"), std::string("z"));
            query2.set_subarray(qc_subarray);
            query2.set_condition(qc.get_value_filter().query_condition());
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

    auto check_empty_query_condition = [&](const CoordinateValueFilter& qc, const std::string& log_note) {
        INFO(log_note);

        // Check valid query condition.
        REQUIRE(qc.is_initialized());

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
        Subarray qc_subarray(*ctx->tiledb_ctx(), array);
        qc_subarray.add_range(0, std::string("a"), std::string("z"));
        query2.set_subarray(qc_subarray);
        query2.set_condition(qc.get_value_filter().query_condition());
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
        CoordinateValueFilter qc(*ctx, dim_names);
        qc.add_column_selection(0, SOMASliceSelection<std::string>("a", "z"));
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("apple"), std::string("fig"));
        check_query_condition(qc, subarray, "Read all values by range.");
    }

    {
        INFO("Invalid points: no values selected - should throw an std::invalid_argument error");
        CoordinateValueFilter qc(*ctx, dim_names);
        std::vector<std::string> points{};
        CHECK_THROWS_AS(qc.add_column_selection(0, SOMAPointSelection<std::string>(points)), std::invalid_argument);
    }

    {
        CoordinateValueFilter qc(*ctx, dim_names);
        qc.add_column_selection(0, SOMASliceSelection<std::string>("ca", "fa"));
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("ca"), std::string("fa"));
        check_query_condition(qc, subarray, "Create by range.");
    }

    {
        CoordinateValueFilter qc(*ctx, dim_names);
        std::vector<std::string> points{"fig", "durian", "banana"};
        qc.add_column_selection(0, SOMAPointSelection<std::string>(points));
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range(0, std::string("banana"), std::string("banana"))
            .add_range(0, std::string("durian"), std::string("durian"))
            .add_range(0, std::string("fig"), std::string("fig"));
        check_query_condition(qc, subarray, "Select by points.");
    }

    {
        CoordinateValueFilter qc(*ctx, dim_names);
        qc.add_column_selection(0, SOMASliceSelection<std::string>("g", "z"));
        check_empty_query_condition(qc, "Create by range - no values selected.");
    }

    {
        CoordinateValueFilter qc(*ctx, dim_names);
        std::vector<std::string> points{"kiwi", "pear", "carrot"};
        qc.add_column_selection(0, SOMAPointSelection<std::string>(points));
        check_empty_query_condition(qc, "Create by points - no values selected.");
    }
}
