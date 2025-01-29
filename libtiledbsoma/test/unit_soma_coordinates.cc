/**
 * @file   unit_soma_coordinates.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for coordinate spaces and coordinate space
 * transformations.
 */

#include "common.h"

TEST_CASE(
    "SOMACoordinateSpace: invalid initialization", "[metadata][spatial]") {
    // Check all three constructors throw for empty coordinate space.
    std::vector<SOMAAxis> empty_axes_{};
    std::vector<std::string> empty_axis_names_{};
    std::vector<std::optional<std::string>> empty_axis_units_{};
    CHECK_THROWS_AS(SOMACoordinateSpace(empty_axes_), TileDBSOMAError);
    CHECK_THROWS_AS(SOMACoordinateSpace(empty_axis_names_), TileDBSOMAError);
    CHECK_THROWS_AS(
        SOMACoordinateSpace(empty_axis_names_, empty_axis_units_),
        TileDBSOMAError);

    // Check all three constructors throw for repeat names.
    std::vector<SOMAAxis> repeat_axes_{
        {"x_axis", "meter"}, {"x_axis", "nanometer"}};
    std::vector<std::string> repeat_axis_names_{"x_axis", "x_axis"};
    std::vector<std::optional<std::string>> repeat_axis_units_{};
    CHECK_THROWS_AS(SOMACoordinateSpace(repeat_axes_), TileDBSOMAError);
    CHECK_THROWS_AS(SOMACoordinateSpace(repeat_axis_names_), TileDBSOMAError);
    CHECK_THROWS_AS(
        SOMACoordinateSpace(repeat_axis_names_, repeat_axis_units_),
        TileDBSOMAError);

    // Check separate axis name/unit constructor throws for mismatched vector
    // sizes.
    std::vector<std::string> axis_names_len3_{"x", "y", "z"};
    std::vector<std::optional<std::string>> axis_units_len2_{"meter", "meter"};
    CHECK_THROWS_AS(
        SOMACoordinateSpace(axis_names_len3_, axis_units_len2_),
        TileDBSOMAError);
    CHECK_THROWS_AS(
        SOMACoordinateSpace(axis_names_len3_, empty_axis_units_),
        TileDBSOMAError);

    // Check valid axis names.
    std::vector<std::string> axis_names_reserved_{"x", "soma_y"};
    std::vector<SOMAAxis> axis_reserved_{{"soma_x", "meter"}, {"y", "meter"}};
    CHECK_THROWS_AS(SOMACoordinateSpace(axis_names_reserved_), TileDBSOMAError);
    CHECK_THROWS_AS(
        SOMACoordinateSpace(axis_names_reserved_, axis_units_len2_),
        TileDBSOMAError);
    CHECK_THROWS_AS(SOMACoordinateSpace(axis_reserved_), TileDBSOMAError);
}

TEST_CASE("SOMACoordinateSpace: compare constructors", "[metadata][spatial]") {
    std::vector<std::string> names1{"x", "y"};
    std::vector<std::optional<std::string>> units1{std::nullopt, std::nullopt};
    SOMACoordinateSpace expected1(
        {SOMAAxis{"x", std::nullopt}, SOMAAxis{"y", std::nullopt}});

    CHECK(SOMACoordinateSpace() == expected1);
    CHECK(SOMACoordinateSpace(names1) == expected1);
    CHECK(SOMACoordinateSpace(names1, units1) == expected1);

    std::vector<std::string> names2{"x1", "x2", "x3"};
    std::vector<std::optional<std::string>> units2{
        "meters", "meters", "meters"};
    SOMACoordinateSpace expected2(
        {SOMAAxis{"x1", "meters"},
         SOMAAxis{"x2", "meters"},
         SOMAAxis{"x3", "meters"}});

    CHECK(SOMACoordinateSpace(names2, units2) == expected2);
}

TEST_CASE("SOMACoordinateSpace: from_metadata", "[metadata][spatial]") {
    // Create the test data.
    std::pair<std::string, SOMACoordinateSpace> test_data = GENERATE(
        std::make_pair(
            "[{\"name\":\"x\",\"unit\":null},{\"name\":\"y\",\"unit\":null}]",
            SOMACoordinateSpace({SOMAAxis{"x"}, SOMAAxis{"y"}})),
        std::make_pair(
            "[{\"name\":\"z\",\"unit\":null},{\"name\":\"t\",\"unit\":"
            "\"hour\"}]",
            SOMACoordinateSpace(
                {SOMAAxis{"z", std::nullopt}, SOMAAxis{"t", "hour"}})),
        std::make_pair(
            "[{\"name\":\"z\",\"unit\":\"meter\"},{\"name\":\"t\",\"unit\":"
            "\"hour\"}]",
            SOMACoordinateSpace(
                {SOMAAxis{"z", "meter"}, SOMAAxis{"t", "hour"}})));
    const auto& metadata = test_data.first;
    const auto& expected_coord_space = test_data.second;

    // Check `from_metadata` creates the expected SOMACoordinateSpace.
    auto actual_coord_space = SOMACoordinateSpace::from_metadata(
        TILEDB_STRING_ASCII,
        static_cast<uint32_t>(metadata.size()),
        static_cast<const void*>(metadata.c_str()));
    REQUIRE(actual_coord_space == expected_coord_space);
}

TEST_CASE("SOMACoordinateSpace: to_string", "[metadata][spatial]") {
    // Create the test data.
    std::pair<std::string, SOMACoordinateSpace> test_data = GENERATE(
        std::make_pair(
            "[{\"name\":\"x\",\"unit\":null},{\"name\":\"y\",\"unit\":null}]",
            SOMACoordinateSpace({SOMAAxis{"x"}, SOMAAxis{"y"}})),
        std::make_pair(
            "[{\"name\":\"z\",\"unit\":null},{\"name\":\"t\",\"unit\":"
            "\"hour\"}]",
            SOMACoordinateSpace(
                {SOMAAxis{"z", std::nullopt}, SOMAAxis{"t", "hour"}})),
        std::make_pair(
            "[{\"name\":\"z\",\"unit\":\"meter\"},{\"name\":\"t\",\"unit\":"
            "\"hour\"}]",
            SOMACoordinateSpace(
                {SOMAAxis{"z", "meter"}, SOMAAxis{"t", "hour"}})));
    const auto& expected_metadata = test_data.first;
    const auto& coord_space = test_data.second;

    auto actual_metadata = coord_space.to_string();
    REQUIRE(actual_metadata == expected_metadata);
}

TEST_CASE(
    "SOMACoordinateSpace: mock metadata round-trip", "[metadata][spatial]") {
    auto coord_space = GENERATE(
        SOMACoordinateSpace(),
        SOMACoordinateSpace(
            {SOMAAxis{"x_axis", "meter"},
             SOMAAxis{"y_axis", "meter"},
             SOMAAxis{"z_axis", "meter"}}),
        SOMACoordinateSpace({SOMAAxis{"z_axis"}, SOMAAxis{"x_axis", "meter"}}));

    auto coord_json = coord_space.to_string();

    auto coord_space_result = SOMACoordinateSpace::from_metadata(
        TILEDB_STRING_ASCII,
        static_cast<uint32_t>(coord_json.size()),
        static_cast<const void*>(coord_json.c_str()));

    REQUIRE(coord_space == coord_space_result);
}
