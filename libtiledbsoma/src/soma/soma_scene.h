/**
 * @file   soma_scene.h
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
 *   This file defines the SOMAScene class.
 */

#ifndef SOMA_SCENE
#define SOMA_SCENE

#include <tiledb/tiledb>

#include "soma_collection.h"

namespace tiledbsoma {

using namespace tiledb;
class SOMAScene : public SOMACollection {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAScene object at the given URI.
     *
     * @param uri URI to create the SOMAScene
     * @param schema TileDB ArraySchema
     * @param platform_config Optional config parameter dictionary
     */
    static void create(
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open a group at the specified URI and return SOMAScene
     * object.
     *
     * @param uri URI of the array
     * @param mode read or write
     * @param ctx TileDB context
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::shared_ptr<SOMAScene> SOMAScene
     */
    static std::unique_ptr<SOMAScene> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    SOMAScene(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMACollection(mode, uri, ctx, timestamp) {
    }

    SOMAScene(const SOMACollection& other)
        : SOMACollection(other) {
    }

    SOMAScene() = delete;
    SOMAScene(const SOMAScene&) = default;
    SOMAScene(SOMAScene&&) = default;
    ~SOMAScene() = default;

    /**
     * @brief Get the collection of imagery data.
     *
     * @return std::shared_ptr<SOMACollection>
     */
    std::shared_ptr<SOMACollection> img();

    /**
     * @brief Get the collection of observation location data.
     *
     * @return std::shared_ptr<SOMACollection>
     */
    std::shared_ptr<SOMACollection> obsl();

    /**
     * @brief Get the collection of collections of variable location data
     * for different measurements.
     *
     * @return std::shared_ptr<SOMACollection>
     */
    std::shared_ptr<SOMACollection> varl();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // A collection of imagery data.
    std::shared_ptr<SOMACollection> img_ = nullptr;

    // A collection of observation location data.
    std::shared_ptr<SOMACollection> obsl_ = nullptr;

    // A collection of collections of variable location data for measurements.
    std::shared_ptr<SOMACollection> varl_ = nullptr;
};
}  // namespace tiledbsoma

#endif  // SOMA_SCENE
