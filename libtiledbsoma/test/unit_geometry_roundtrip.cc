/** * @file   unit_geometry_roundtrip.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for testing the geometry utilities
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
#include "../src/geometry/operators/io/read.h"
#include "../src/geometry/operators/io/write.h"

using namespace tiledbsoma;
using namespace tiledbsoma::geometry;
using namespace Catch::Matchers;

#ifndef TILEDBSOMA_SOURCE_ROOT
#define TILEDBSOMA_SOURCE_ROOT "not_defined"
#endif

TEST_CASE("Geometry: Point roundtrip") {
    Point point(GENERATE(1.0, 5.2), GENERATE(3.0, 1.3));

    BinaryBuffer buffer = to_wkb(GenericGeometry(point));
    CHECK(buffer.size() == 21);

    GenericGeometry parsedGeometry = from_wkb<GenericGeometry>(buffer);

    CHECK(std::holds_alternative<Point>(parsedGeometry));
    Point parsedPoint = std::get<Point>(parsedGeometry);
    CHECK(point.x == parsedPoint.x);
    CHECK(point.y == parsedPoint.y);
}

TEST_CASE("Geometry: LineString roundtrip") {
    LineString linestring(std::vector<BasePoint>(
        {BasePoint(5.2, 1.3), BasePoint(1.0, 3.0), BasePoint(1.0, 1.4)}));

    BinaryBuffer buffer = to_wkb(GenericGeometry(linestring));
    CHECK(buffer.size() == 57);

    GenericGeometry parsedGeometry = from_wkb<GenericGeometry>(buffer);

    CHECK(std::holds_alternative<LineString>(parsedGeometry));
    LineString parsedLineString = std::get<LineString>(parsedGeometry);

    CHECK(parsedLineString.points.size() == 3);
    for (uint32_t i = 0; i < linestring.points.size(); ++i) {
        CHECK(linestring.points[i].x == parsedLineString.points[i].x);
        CHECK(linestring.points[i].y == parsedLineString.points[i].y);
    }
}

TEST_CASE("Geometry: Polygon roundtrip") {
    Polygon polygon(
        std::vector<BasePoint>(
            {BasePoint(5.2, 1.3), BasePoint(1.0, 3.0), BasePoint(1.0, 1.4)}),
        std::vector<std::vector<BasePoint>>({std::vector<BasePoint>(
            {BasePoint(5.234, 1.3), BasePoint(13.0, 3.0)})}));

    BinaryBuffer buffer = to_wkb(GenericGeometry(polygon));
    CHECK(buffer.size() == 97);

    GenericGeometry parsedGeometry = from_wkb<GenericGeometry>(buffer);

    CHECK(std::holds_alternative<Polygon>(parsedGeometry));

    Polygon parsedPolygon = std::get<Polygon>(parsedGeometry);

    CHECK(parsedPolygon.exteriorRing.size() == 3);
    CHECK(parsedPolygon.interiorRings.size() == 1);
    CHECK(parsedPolygon.interiorRings[0].size() == 2);

    for (uint32_t i = 0; i < polygon.exteriorRing.size(); ++i) {
        CHECK(parsedPolygon.exteriorRing[i].x == polygon.exteriorRing[i].x);
        CHECK(parsedPolygon.exteriorRing[i].y == polygon.exteriorRing[i].y);
    }

    for (uint32_t i = 0; i < polygon.interiorRings.size(); ++i) {
        for (uint32_t j = 0; j < polygon.interiorRings[i].size(); ++j) {
            CHECK(
                parsedPolygon.interiorRings[i][j].x ==
                polygon.interiorRings[i][j].x);
            CHECK(
                parsedPolygon.interiorRings[i][j].y ==
                polygon.interiorRings[i][j].y);
        }
    }
}

TEST_CASE("Geometry: Multipoint roundtrip") {
    MultiPoint multi_point(std::vector<Point>(
        {Point(GENERATE(1.0, 5.2), 1.3), Point(3.14, GENERATE(3.0, 0))}));

    BinaryBuffer buffer = to_wkb(GenericGeometry(multi_point));
    CHECK(buffer.size() == 51);

    GenericGeometry parsedGeometry = from_wkb<GenericGeometry>(buffer);

    CHECK(std::holds_alternative<MultiPoint>(parsedGeometry));

    MultiPoint parsed_multi_point = std::get<MultiPoint>(parsedGeometry);

    CHECK(parsed_multi_point.points.size() == multi_point.points.size());

    for (uint32_t i = 0; i < multi_point.points.size(); ++i) {
        CHECK(multi_point.points[i].x == parsed_multi_point.points[i].x);
        CHECK(multi_point.points[i].y == parsed_multi_point.points[i].y);
    }
}

TEST_CASE("Geometry: MultiLineString roundtrip") {
    MultiLineString multi_linestring(std::vector<LineString>(
        {LineString(std::vector<BasePoint>(
             {BasePoint(5.2, 1.3), BasePoint(1.0, 3.0), BasePoint(1.0, 1.4)})),
         LineString(std::vector<BasePoint>(
             {BasePoint(0, 100), BasePoint(3.15, 1.67)}))}));

    BinaryBuffer buffer = to_wkb(GenericGeometry(multi_linestring));
    CHECK(buffer.size() == 107);

    GenericGeometry parsedGeometry = from_wkb<GenericGeometry>(buffer);

    CHECK(std::holds_alternative<MultiLineString>(parsedGeometry));

    MultiLineString parsed_multi_linestring = std::get<MultiLineString>(
        parsedGeometry);

    CHECK(
        parsed_multi_linestring.linestrings.size() ==
        multi_linestring.linestrings.size());

    for (uint32_t i = 0; i < multi_linestring.linestrings.size(); ++i) {
        CHECK(
            parsed_multi_linestring.linestrings[i].points.size() ==
            multi_linestring.linestrings[i].points.size());
        for (uint32_t j = 0; j < multi_linestring.linestrings[i].points.size();
             ++j) {
            CHECK(
                multi_linestring.linestrings[i].points[j].x ==
                parsed_multi_linestring.linestrings[i].points[j].x);
            CHECK(
                multi_linestring.linestrings[i].points[j].y ==
                parsed_multi_linestring.linestrings[i].points[j].y);
        }
    }
}

TEST_CASE("Geometry: MultiPolygon roundtrip") {
    MultiPolygon multi_polygon(std::vector<Polygon>(
        {Polygon(
             std::vector<BasePoint>(
                 {BasePoint(5.2, 1.3),
                  BasePoint(1.0, 3.0),
                  BasePoint(1.0, 1.4)}),
             std::vector<std::vector<BasePoint>>({std::vector<BasePoint>(
                 {BasePoint(5.234, 1.3), BasePoint(13.0, 3.0)})})),
         Polygon(std::vector<BasePoint>(
             {BasePoint(5.2, 6.3),
              BasePoint(1.3664, 3.0),
              BasePoint(1, 10)}))}));

    BinaryBuffer buffer = to_wkb(GenericGeometry(multi_polygon));
    CHECK(buffer.size() == 167);

    GenericGeometry parsedGeometry = from_wkb<GenericGeometry>(buffer);

    CHECK(std::holds_alternative<MultiPolygon>(parsedGeometry));

    MultiPolygon parsed_multi_polygon = std::get<MultiPolygon>(parsedGeometry);

    CHECK(
        parsed_multi_polygon.polygons.size() == multi_polygon.polygons.size());

    for (uint32_t i = 0; i < multi_polygon.polygons.size(); ++i) {
        CHECK(
            multi_polygon.polygons[i].exteriorRing.size() ==
            parsed_multi_polygon.polygons[i].exteriorRing.size());

        for (uint32_t j = 0; j < multi_polygon.polygons[i].exteriorRing.size();
             ++j) {
            CHECK(
                multi_polygon.polygons[i].exteriorRing[j].x ==
                parsed_multi_polygon.polygons[i].exteriorRing[j].x);
            CHECK(
                multi_polygon.polygons[i].exteriorRing[j].y ==
                parsed_multi_polygon.polygons[i].exteriorRing[j].y);
        }

        CHECK(
            multi_polygon.polygons[i].interiorRings.size() ==
            parsed_multi_polygon.polygons[i].interiorRings.size());

        for (uint32_t j = 0; j < multi_polygon.polygons[i].interiorRings.size();
             ++j) {
            CHECK(
                multi_polygon.polygons[i].interiorRings[j].size() ==
                parsed_multi_polygon.polygons[i].interiorRings[j].size());

            for (uint32_t k = 0;
                 k < multi_polygon.polygons[i].interiorRings[j].size();
                 ++k) {
                CHECK(
                    multi_polygon.polygons[i].interiorRings[j][k].x ==
                    parsed_multi_polygon.polygons[i].interiorRings[j][k].x);
                CHECK(
                    multi_polygon.polygons[i].interiorRings[j][k].y ==
                    parsed_multi_polygon.polygons[i].interiorRings[j][k].y);
            }
        }
    }
}
