/**
 * @file   unit_soma_scene.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
