/**
 * @file  test_soma_column_selection.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages tests for the SOMA column selection classes.
 */

#include "common.h"

TEMPLATE_TEST_CASE(
    "SliceSelection: numeric (signed)",
    "[SOMAColumnSelection][SliceSelection]",
    int8_t,
    int16_t,
    int32_t,
    int64_t,
    float_t,
    double_t) {
    SliceSelection<TestType> slice1(0, 10);
    CHECK(slice1.start == 0);
    CHECK(slice1.stop == 10);
    CHECK(slice1.has_overlap({-100, 100}));
    CHECK(slice1.has_overlap({2, 3}));
    CHECK(slice1.has_overlap({-10, 0}));
    CHECK(slice1.has_overlap({10, 11}));
    CHECK(!slice1.has_overlap({-10, -1}));
    CHECK(!slice1.has_overlap({11, 300}));

    SliceSelection<TestType> slice2(std::make_pair(-10, 10));
    CHECK(slice2.start == -10);
    CHECK(slice2.stop == 10);
    CHECK(slice2.has_overlap({-100, 100}));
    CHECK(slice2.has_overlap({-4, 4}));
    CHECK(slice2.has_overlap({-100, -10}));
    CHECK(slice2.has_overlap({10, 10}));
    CHECK(!slice2.has_overlap({-100, -11}));
    CHECK(!slice2.has_overlap({11, 100}));

    SliceSelection<TestType> slice3(std::make_pair(8, 8));
    CHECK(slice3.start == 8);
    CHECK(slice3.stop == 8);
    CHECK(slice3.has_overlap({-100, 100}));
    CHECK(slice3.has_overlap({0, 8}));
    CHECK(slice3.has_overlap({8, 10}));
    CHECK(slice3.has_overlap({8, 8}));
    CHECK(!slice3.has_overlap({-100, -11}));
    CHECK(!slice3.has_overlap({11, 100}));

    SliceSelection<TestType> slice4(std::nullopt, std::nullopt);
    CHECK(!slice4.start.has_value());
    CHECK(!slice4.stop.has_value());
    CHECK(slice4.has_overlap({-100, 100}));

    SliceSelection<TestType> slice5(10, std::nullopt);
    CHECK(slice5.start == 10);
    CHECK(!slice5.stop.has_value());
    CHECK(slice5.has_overlap({-100, 100}));
    CHECK(slice5.has_overlap({50, 100}));
    CHECK(!slice5.has_overlap({-100, 0}));

    SliceSelection<TestType> slice6(std::nullopt, 10);
    CHECK(!slice6.start.has_value());
    CHECK(slice6.stop == 10);
    CHECK(slice6.has_overlap({-100, 100}));
    CHECK(slice6.has_overlap({-100, 0}));
    CHECK(!slice6.has_overlap({50, 100}));

    CHECK_THROWS_AS(SliceSelection<TestType>(10, -10), std::invalid_argument);
}

TEMPLATE_TEST_CASE(
    "SliceSelection: numeric (unsigned)",
    "[SOMAColumnSelection][SliceSelection]",
    uint8_t,
    uint16_t,
    uint32_t,
    uint64_t) {
    SliceSelection<TestType> slice1(0, 10);
    CHECK(slice1.start == 0);
    CHECK(slice1.stop == 10);
    CHECK(slice1.has_overlap({0, 100}));
    CHECK(slice1.has_overlap({2, 3}));
    CHECK(slice1.has_overlap({0, 0}));
    CHECK(slice1.has_overlap({10, 11}));
    CHECK(!slice1.has_overlap({11, 300}));

    SliceSelection<TestType> slice2(std::make_pair(4, 10));
    CHECK(slice2.start == 4);
    CHECK(slice2.stop == 10);
    CHECK(slice2.has_overlap({0, 100}));
    CHECK(slice2.has_overlap({0, 4}));
    CHECK(slice2.has_overlap({10, 10}));
    CHECK(!slice2.has_overlap({0, 3}));
    CHECK(!slice2.has_overlap({11, 100}));

    SliceSelection<TestType> slice3(std::make_pair(8, 8));
    CHECK(slice3.start == 8);
    CHECK(slice3.stop == 8);
    CHECK(slice3.has_overlap({0, 8}));
    CHECK(slice3.has_overlap({8, 10}));
    CHECK(slice3.has_overlap({8, 8}));
    CHECK(slice3.has_overlap({4, 10}));
    CHECK(!slice3.has_overlap({0, 7}));
    CHECK(!slice3.has_overlap({9, 100}));

    SliceSelection<TestType> slice4(std::nullopt, std::nullopt);
    CHECK(!slice4.start.has_value());
    CHECK(!slice4.stop.has_value());
    CHECK(slice4.has_overlap({0, 100}));

    SliceSelection<TestType> slice5(10, std::nullopt);
    CHECK(slice5.start == 10);
    CHECK(!slice5.stop.has_value());
    CHECK(slice5.has_overlap({0, 100}));
    CHECK(slice5.has_overlap({50, 100}));
    CHECK(!slice5.has_overlap({0, 8}));

    SliceSelection<TestType> slice6(std::nullopt, 10);
    CHECK(!slice6.start.has_value());
    CHECK(slice6.stop == 10);
    CHECK(slice6.has_overlap({0, 100}));
    CHECK(slice6.has_overlap({0, 8}));
    CHECK(!slice6.has_overlap({50, 100}));

    CHECK_THROWS_AS(SliceSelection<TestType>(10, 0), std::invalid_argument);
}

TEMPLATE_TEST_CASE(
    "PointSelection: numeric (signed)",
    "[SOMAColumnSelection][PointSelection]",
    int8_t,
    int16_t,
    int32_t,
    int64_t,
    float_t,
    double_t) {
    std::vector<TestType> values{1, 54, -3, 17};
    PointSelection<TestType> points1(values);
    CHECK(points1.is_subset({-3, 54}));
    CHECK(points1.is_subset({-3, 54}));
    CHECK(!points1.is_subset({1, 17}));
    CHECK(!points1.is_subset({-5, 40}));
    CHECK(!points1.is_subset({-2, 56}));

    PointSelection<TestType> points2(std::span<TestType>{});
    CHECK(points2.is_subset({-100, 100}));
}

TEMPLATE_TEST_CASE(
    "PointSelection: numeric (unsigned)",
    "[SOMAColumnSelection][PointSelection]",
    uint8_t,
    uint16_t,
    uint32_t,
    uint64_t) {
    std::vector<TestType> points{1, 54, 2, 17};
    PointSelection<TestType> points1(points);
    CHECK(points1.is_subset({1, 54}));
    CHECK(points1.is_subset({0, 100}));
    CHECK(!points1.is_subset({1, 17}));
    CHECK(!points1.is_subset({0, 40}));
    CHECK(!points1.is_subset({2, 56}));

    PointSelection<TestType> points2(std::span<TestType>{});
    CHECK(points2.is_subset({-100, 100}));
}
