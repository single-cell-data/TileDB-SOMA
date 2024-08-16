/**
 * @file   soma_image2d_collection.h
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

#ifndef SOMA_IMAGE2D_COLLECTION
#define SOMA_IMAGE2D_COLLECTION

#include <tiledb/tiledb>

#include "soma_collection.h"

namespace tiledbsoma {

using namespace tiledb;
class SOMAImage2DCollection : public SOMACollection {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAImage2DCollection object at the given URI.
     *
     * @param uri URI to create the SOMAImage2DCollection
     * @param schema TileDB ArraySchema
     * @param platform_config Optional config parameter dictionary
     */
    static void create(
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open a group at the specified URI and return SOMAImage2DCollection
     * object.
     *
     * @param uri URI of the array
     * @param mode read or write
     * @param ctx TileDB context
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::shared_ptr<SOMAImage2DCollection> SOMAImage2DCollection
     */
    static std::unique_ptr<SOMAImage2DCollection> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    SOMAImage2DCollection(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMACollection(mode, uri, ctx, timestamp) {
    }

    SOMAImage2DCollection(const SOMACollection& other)
        : SOMACollection(other) {
    }

    SOMAImage2DCollection() = delete;
    SOMAImage2DCollection(const SOMAImage2DCollection&) = default;
    SOMAImage2DCollection(SOMAImage2DCollection&&) = default;
    ~SOMAImage2DCollection() = default;

   private:
    //===================================================================
    //= private non-static
    //===================================================================
};
}  // namespace tiledbsoma

#endif  // SOMA_IMAGE2D_COLLECTION
