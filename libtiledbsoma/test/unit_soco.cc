/**
 * @file   unit_soco.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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
 * This file manages unit tests for soma collection objects
 */

#include <catch2/catch_test_macros.hpp>
#include <tiledb/tiledb>
#include <tiledbsoma/tiledbsoma>

#ifndef TILEDBSOMA_SOURCE_ROOT
#define TILEDBSOMA_SOURCE_ROOT "not_defined"
#endif

static const std::string root = TILEDBSOMA_SOURCE_ROOT;
static const std::string soco_uri = root + "/test/soco";

using namespace tiledb;
using namespace tiledbsoma;

TEST_CASE("SOCO: Open arrays") {
    Config config;
    // config.logging_level"] = "5";

    auto soco = SOMACollection::open(soco_uri, config);
    auto soma_uris = soco->list_somas();
    REQUIRE(soma_uris.size() == 2);

    for (const auto& [name, uri] : soma_uris) {
        (void)name;
        auto soma = SOMA::open(uri);
    }
}
