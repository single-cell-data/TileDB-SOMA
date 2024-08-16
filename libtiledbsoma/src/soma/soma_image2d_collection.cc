/**
 * @file   soma_image2d_collection.cc
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
 *   This file defines the SOMAImage2DCollection class.
 */

#include "soma_image2d_collection.h"
#include "soma_collection.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAImage2DCollection::create(
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    try {
        std::filesystem::path image_uri(uri);
        SOMAGroup::create(
            ctx, image_uri.string(), "SOMAImage2DCollection", timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAImage2DCollection> SOMAImage2DCollection::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    try {
        return std::make_unique<SOMAImage2DCollection>(
            mode, uri, ctx, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

}  // namespace tiledbsoma
