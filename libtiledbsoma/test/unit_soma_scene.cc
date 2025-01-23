/**
 * @file   unit_soma_scene.cc
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
 * This file manages unit tests for the SOMAScene class
 */

#include "common.h"

TEST_CASE("SOMAScene: basic", "[scene][spatial]") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri{"mem://unit-test-scene-basic"};

    SOMAScene::create(uri, ctx, std::nullopt);
    auto soma_scene = SOMAScene::open(uri, OpenMode::read, ctx, std::nullopt);
    CHECK(soma_scene->uri() == uri);
    CHECK(soma_scene->ctx() == ctx);
    CHECK(soma_scene->type() == "SOMAScene");
    CHECK(soma_scene->has_metadata("soma_encoding_version"));
    CHECK(soma_scene->has_metadata("soma_spatial_encoding_version"));
    CHECK(not soma_scene->coordinate_space().has_value());
    soma_scene->close();
}

TEST_CASE("SOMAScene: with coordinates", "[scene][spatial]") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri{"mem://unit-test-scene-coords"};

    SOMACoordinateSpace coord_space{};

    SOMAScene::create(uri, ctx, coord_space, std::nullopt);

    auto soma_scene = SOMAScene::open(uri, OpenMode::read, ctx, std::nullopt);
    CHECK(soma_scene->uri() == uri);
    CHECK(soma_scene->ctx() == ctx);
    CHECK(soma_scene->type() == "SOMAScene");
    CHECK(soma_scene->has_metadata("soma_encoding_version"));
    CHECK(soma_scene->has_metadata("soma_spatial_encoding_version"));
    auto scene_coord_space = soma_scene->coordinate_space();
    REQUIRE(scene_coord_space.has_value());

    CHECK(scene_coord_space.value() == coord_space);
}
