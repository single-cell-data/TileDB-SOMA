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
 * Test that checks SOMAQueryCondition is grabbing the correct values.
 *
 * This is an ideal test to parametrize with rapidtest when we add it to SOMA.
 */
TEST_CASE("Test SOMAQueryCondition on SparseArray", "[SOMAQueryCondition][SOMASparseNDArray]") {
    // Create a TileDB sparse array with 3 integer dimensions.
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
    auto check_query_condition =
        [&](const SOMAQueryCondition& qc, const Subarray& subarray, const std::string& log_note) {
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

    auto check_empty_query_condition = [&](const SOMAQueryCondition& qc, const std::string& log_note) {
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
        SOMAQueryCondition qc{};
        CHECK(!qc.is_initialized());
    }

    // Full region by range.
    {
        auto qc = SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 0, 15);
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 0, 15);
        check_query_condition(qc, subarray, "Read all values by range.");
    }

    // Empty region: invalid range.
    {
        auto qc = SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 3, 2);
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        check_empty_query_condition(qc, "Invalid range: expect no values.");
    }

    // Empty region: out-of-bounds range.
    {
        auto qc = SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 25, 27);
        CHECK(qc.is_initialized());
        check_empty_query_condition(qc, "Range out of bounds: expect no values.");
    }

    // Empty region: out-of-bounds points.
    {
        auto qc = SOMAQueryCondition::create_from_points<int64_t>(*tiledb_ctx, "soma_dim_0", {25, 27, 19, 30});
        CHECK(qc.is_initialized());
        check_empty_query_condition(qc, "Range out of bounds: expect no values.");
    }

    // Region [1:2] by ranges.
    {
        auto qc = SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 1, 2);
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 1, 2);
        check_query_condition(qc, subarray, "Select by range on dim 0.");
    }

    // Region [0,11,13] by points.
    {
        auto qc = SOMAQueryCondition::create_from_points<int64_t>(*tiledb_ctx, "soma_dim_0", {0, 1, 13});
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 0, 1).add_range<int64_t>("soma_dim_0", 13, 13);
        check_query_condition(qc, subarray, "Select by points on dim 0.");
    }

    // Test "combine_with_and".
    {
        INFO("Check `combine_with_and` where neither value is initialized.");
        auto qc = SOMAQueryCondition().combine_with_and(SOMAQueryCondition());
        CHECK(!qc.is_initialized());
    }
    {
        auto qc = SOMAQueryCondition().combine_with_and(
            SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 0, 0));
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 0, 0);
        check_query_condition(qc, subarray, "Check 'combine_with_and' second value initialized.");
    }
    {
        auto qc = SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 0, 2)
                      .combine_with_and(SOMAQueryCondition());
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 0, 2);
        check_query_condition(qc, subarray, "Check 'combine_with_and' first value initialized.");
    }
    {
        auto qc = SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 1, 3)
                      .combine_with_and(
                          SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 0, 2));
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 1, 2);
        check_query_condition(qc, subarray, "Check 'combine_with_and' both values initialized.");
    }

    // Test "combine with or".
    {
        INFO("Check `combine_with_and` where neither value is initialized.");
        auto qc = SOMAQueryCondition().combine_with_or(SOMAQueryCondition());
        CHECK(!qc.is_initialized());
    }
    {
        auto qc = SOMAQueryCondition().combine_with_or(
            SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 0, 2));
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 0, 2);
        check_query_condition(qc, subarray, "Check 'combine_with_or' second value initialized.");
    }
    {
        auto qc = SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 0, 2)
                      .combine_with_or(SOMAQueryCondition());
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 0, 2);
        check_query_condition(qc, subarray, "Check 'combine_with_or' first value initialized.");
    }
    {
        auto qc = SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 1, 3)
                      .combine_with_or(
                          SOMAQueryCondition::create_from_range<int64_t>(*tiledb_ctx, "soma_dim_0", 11, 13));
        CHECK(qc.is_initialized());
        Subarray subarray(*tiledb_ctx, array);
        subarray.add_range<int64_t>("soma_dim_0", 1, 3).add_range<int64_t>("soma_dim_0", 11, 13);
        check_query_condition(qc, subarray, "Check 'combine_with_or' both values initialized.");
    }
}

/**
 * Test that checks SOMArdQueryCondition is grabbing the correct values.
 *
 * This is an ideal test to parametrize with rapidtest when we add it to SOMA.
 */
TEST_CASE("Test SOMACoordQueryCondition on SparseArray", "[SOMACoordQueryCondition][SOMASparseNDArray]") {
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
        [&](const SOMACoordQueryCondition& qc, const Subarray& subarray, const std::string& log_note) {
            INFO(log_note);

            // Create a query for the entire array.
            std::vector<int64_t> expected_coords_0(16);
            std::vector<int64_t> expected_coords_1(16);
            std::vector<int32_t> expected_data(16);

            // Create a query for the entire array.
            std::vector<int64_t> actual_coords_0(16);
            std::vector<int64_t> actual_coords_1(16);
            std::vector<int32_t> actual_data(16);

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
            qc_subarray.add_range<int64_t>(0, 0, 3).add_range<int64_t>(1, 0, 3);
            query2.set_subarray(qc_subarray);
            query2.set_condition(qc.query_condition());
            query2.submit();

            // Check results.
            auto expected_result_num = static_cast<int64_t>(query1.result_buffer_elements()["soma_data"].second);
            auto actual_result_num = static_cast<int64_t>(query2.result_buffer_elements()["soma_data"].second);
            CHECK(expected_result_num == actual_result_num);
            CHECK_THAT(actual_coords_0, Catch::Matchers::Equals(expected_coords_0));
            CHECK_THAT(actual_coords_1, Catch::Matchers::Equals(expected_coords_1));
            CHECK_THAT(actual_data, Catch::Matchers::Equals(expected_data));
        };

    auto check_empty_query_condition = [&](const SOMACoordQueryCondition& qc, const std::string& log_note) {
        INFO(log_note);

        // Create a query for the entire array.
        std::vector<int64_t> expected_coords_0(16);
        std::vector<int64_t> expected_coords_1(16);
        std::vector<int32_t> expected_data(16);

        // Create a query for the entire array.
        std::vector<int64_t> actual_coords_0(16);
        std::vector<int64_t> actual_coords_1(16);
        std::vector<int32_t> actual_data(16);

        // Query with query condition.
        Query query2(*ctx->tiledb_ctx(), array);
        query2.set_layout(TILEDB_ROW_MAJOR);
        query2.set_data_buffer("soma_dim_0", actual_coords_0);
        query2.set_data_buffer("soma_dim_1", actual_coords_1);
        query2.set_data_buffer("soma_data", actual_data);
        Subarray qc_subarray(*ctx->tiledb_ctx(), array);
        qc_subarray.add_range<int64_t>(0, 0, 3).add_range<int64_t>(1, 0, 3);
        query2.set_subarray(qc_subarray);
        query2.set_condition(qc.query_condition());
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
    {
        SOMACoordQueryCondition qc(*ctx);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_range<int64_t>("soma_dim_0", 0, 3).add_range<int64_t>("soma_dim_1", 0, 2);
        subarray.add_range<int64_t>("soma_dim_0", 0, 3).add_range<int64_t>("soma_dim_1", 0, 2);
        check_query_condition(qc, subarray, "Read all values by range.");
    }

    // Empty region: invalid range.
    {
        SOMACoordQueryCondition qc(*ctx);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_range<int64_t>("soma_dim_0", 3, 2);
        check_empty_query_condition(qc, "Invalid range: expect no values.");
    }

    // Empty region: out-of-bounds range.
    {
        SOMACoordQueryCondition qc(*ctx);
        qc.add_range<int64_t>("soma_dim_0", 5, 7);
        check_empty_query_condition(qc, "Range out of bounds: expect no values.");
    }

    // Empty region: out-of-bounds points.
    {
        SOMACoordQueryCondition qc(*ctx);
        qc.add_points<int64_t>("soma_dim_0", {5, 7, 11, 10});
        check_empty_query_condition(qc, "Range out of bounds: expect no values.");
    }

    // Region [1:2]x[:] by ranges.
    {
        SOMACoordQueryCondition qc(*ctx);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_range<int64_t>("soma_dim_0", 1, 2);
        subarray.add_range<int64_t>("soma_dim_0", 1, 2);
        check_query_condition(qc, subarray, "Select by range on dim 0.");
    }

    // Region [:]x[1:2] by ranges.
    {
        SOMACoordQueryCondition qc(*ctx);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_range<int64_t>("soma_dim_1", 1, 2);
        subarray.add_range<int64_t>("soma_dim_1", 1, 2);
        check_query_condition(qc, subarray, "Select by range on dim 1.");
    }

    // Region [0,1,3]x[:] by points (ordered).
    {
        SOMACoordQueryCondition qc(*ctx);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_points<int64_t>("soma_dim_0", {0, 1, 3});
        subarray.add_range<int64_t>("soma_dim_0", 0, 1).add_range<int64_t>("soma_dim_0", 3, 3);
        check_query_condition(qc, subarray, "Select by points on dim 0 (ordered).");
    }

    // Region [0,1,3]x[:] by points (unordered).
    {
        SOMACoordQueryCondition qc(*ctx);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_points<int64_t>("soma_dim_0", {3, 0, 1});
        subarray.add_range<int64_t>("soma_dim_0", 0, 1).add_range<int64_t>("soma_dim_0", 3, 3);
        check_query_condition(qc, subarray, "Select by points on dim 0 (unordered).");
    }

    // Region [0,1,3]x[:] by points (multiple conditions).
    {
        SOMACoordQueryCondition qc(*ctx);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_points<int64_t>("soma_dim_0", {3})
            .add_points<int64_t>("soma_dim_0", {0})
            .add_points<int64_t>("soma_dim_0", {1});
        subarray.add_range<int64_t>("soma_dim_0", 0, 1).add_range<int64_t>("soma_dim_0", 3, 3);
        check_query_condition(qc, subarray, "Select by points on dim 0 (multiple conditions).");
    }

    // Region [:]x[0,2] by points.
    {
        SOMACoordQueryCondition qc(*ctx);
        Subarray subarray(*ctx->tiledb_ctx(), array);
        qc.add_points<int64_t>("soma_dim_1", {0, 2});
        subarray.add_range<int64_t>("soma_dim_1", 0, 0).add_range<int64_t>("soma_dim_1", 2, 2);
        check_query_condition(qc, subarray, "Select by points on dim 1.");
    }
}
