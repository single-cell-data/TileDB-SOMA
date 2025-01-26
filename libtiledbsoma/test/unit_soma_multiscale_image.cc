/**
 * @file   unit_soma_multiscale_image.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMAMultiscaleImage class
 */
#include "common.h"

TEST_CASE("SOMAMultiscaleImage: basic", "[multiscale_image][spatial]") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-multiscale-image-basic";

    SOMACoordinateSpace coord_space{};
    SOMAMultiscaleImage::create(uri, ctx, coord_space, std::nullopt);
    auto soma_image = SOMAMultiscaleImage::open(
        uri, OpenMode::read, ctx, std::nullopt);
    REQUIRE(soma_image->uri() == uri);
    REQUIRE(soma_image->ctx() == ctx);
    REQUIRE(soma_image->type() == "SOMAMultiscaleImage");
    REQUIRE(soma_image->coordinate_space() == coord_space);
    soma_image->close();
}
