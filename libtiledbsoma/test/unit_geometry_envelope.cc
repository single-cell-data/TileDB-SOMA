/**
 * @file   unit_geometry_envelope.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for testing the geometry envelope computation
 */

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_predicate.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <utility>

#include "../src/geometry/geometry.h"
#include "../src/geometry/operators/envelope.h"

using namespace tiledbsoma;
using namespace tiledbsoma::geometry;
using namespace Catch::Matchers;

#ifndef TILEDBSOMA_SOURCE_ROOT
#define TILEDBSOMA_SOURCE_ROOT "not_defined"
#endif

TEST_CASE("Geometry: Envelope") {
    Point point(10, 20);
    LineString linestring(
        std::vector<BasePoint>({BasePoint(2, 2), BasePoint(3, 4)}));
    Polygon polygon(std::vector<BasePoint>(
        {BasePoint(0, 0), BasePoint(1, 0), BasePoint(0, 1)}));
    MultiPoint multi_point(std::vector<Point>({Point(0, 0), Point(1, 1)}));
    MultiLineString multi_linestring(std::vector<LineString>(
        {LineString(std::vector<BasePoint>({BasePoint(2, 2), BasePoint(3, 4)})),
         LineString(
             std::vector<BasePoint>({BasePoint(0, -2), BasePoint(3, 0)}))}));
    MultiPolygon multi_polygon(std::vector<Polygon>(
        {Polygon(std::vector<BasePoint>(
             {BasePoint(0, 1), BasePoint(1, 1), BasePoint(-10, 10)})),
         Polygon(std::vector<BasePoint>(
             {BasePoint(-2, -2),
              BasePoint(-2, 2),
              BasePoint(2, 2),
              BasePoint(2, -2)}))}));
    GeometryCollection collection(
        {point,
         linestring,
         polygon,
         multi_point,
         multi_linestring,
         multi_polygon});

    Envelope point_envelope = envelope(point);
    CHECK(point_envelope.range.at(0).first == point.x);
    CHECK(point_envelope.range.at(0).second == point.x);
    CHECK(point_envelope.range.at(1).first == point.y);
    CHECK(point_envelope.range.at(1).second == point.y);

    Envelope linestring_envelope = envelope(linestring);
    CHECK(linestring_envelope.range.at(0).first == 2);
    CHECK(linestring_envelope.range.at(0).second == 3);
    CHECK(linestring_envelope.range.at(1).first == 2);
    CHECK(linestring_envelope.range.at(1).second == 4);

    Envelope polygon_envelope = envelope(polygon);
    CHECK(polygon_envelope.range.at(0).first == 0);
    CHECK(polygon_envelope.range.at(0).second == 1);
    CHECK(polygon_envelope.range.at(1).first == 0);
    CHECK(polygon_envelope.range.at(1).second == 1);

    Envelope multi_point_envelope = envelope(multi_point);
    CHECK(multi_point_envelope.range.at(0).first == 0);
    CHECK(multi_point_envelope.range.at(0).second == 1);
    CHECK(multi_point_envelope.range.at(1).first == 0);
    CHECK(multi_point_envelope.range.at(1).second == 1);

    Envelope multi_linestring_envelope = envelope(multi_linestring);
    CHECK(multi_linestring_envelope.range.at(0).first == 0);
    CHECK(multi_linestring_envelope.range.at(0).second == 3);
    CHECK(multi_linestring_envelope.range.at(1).first == -2);
    CHECK(multi_linestring_envelope.range.at(1).second == 4);

    Envelope multi_polygon_envelope = envelope(multi_polygon);
    CHECK(multi_polygon_envelope.range.at(0).first == -10);
    CHECK(multi_polygon_envelope.range.at(0).second == 2);
    CHECK(multi_polygon_envelope.range.at(1).first == -2);
    CHECK(multi_polygon_envelope.range.at(1).second == 10);

    Envelope collection_envelope = envelope(collection);
    CHECK(collection_envelope.range.at(0).first == -10);
    CHECK(collection_envelope.range.at(0).second == 10);
    CHECK(collection_envelope.range.at(1).first == -2);
    CHECK(collection_envelope.range.at(1).second == 20);
}
