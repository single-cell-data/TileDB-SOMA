/**
 * @file   soma_collection.h
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
 *   This declares the soma collection classes
 */

#ifndef SOMA_COLLECTION_H
#define SOMA_COLLECTION_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>
#include <tiledbsoma/tiledbsoma>

#include "tiledbsoma/soma_collection_query.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMACollection {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open a SOMACollection at the specified URI and return a
     * SOMACollection object.
     *
     * @param uri URI of SOMACollection
     * @param uri TileDB context
     * @return SOMACollection object
     */
    static std::unique_ptr<SOMACollection> open(
        std::string_view uri,
        std::shared_ptr<Context> ctx = std::make_shared<Context>());

    /**
     * @brief Open a SOMACollection at the specified URI and return a
     * SOMACollection object.
     *
     * @param uri URI of SOMACollection
     * @param config TileDB config
     * @return SOMACollection object
     */
    static std::unique_ptr<SOMACollection> open(
        std::string_view uri, const Config& config);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMACollection object
     *
     * @param uri URI of SOMACollection
     */
    SOMACollection(std::string_view uri, std::shared_ptr<Context> ctx);

    /**
     * @brief Return a map of hierarchical SOMA names to SOMA URIs for all
     * SOMAs in the SOMACollection. This includes nested SOMACollections.
     *
     * @return std::unordered_map<std::string, std::string> Map of SOMA name to
     * SOMA URI
     */
    std::unordered_map<std::string, std::string> list_somas();

    std::unique_ptr<SOMACollectionQuery> query() {
        return std::make_unique<SOMACollectionQuery>(this);
    }

    std::unordered_map<std::string, std::shared_ptr<SOMA>> get_somas();

    /**
     * @brief Return TileDB context.
     *
     * @return std::shared_ptr<Context> Context.
     */
    std::shared_ptr<Context> context() {
        return ctx_;
    }

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMACollection URI
    std::string uri_;

    // Map of hierarchical SOMA name to SOMA URI
    std::unordered_map<std::string, std::string> soma_uri_map_;

    // Map of hierarchical SOMA name to SOMA
    std::unordered_map<std::string, std::shared_ptr<SOMA>> soma_map_;

    /**
     * @brief Walk the TileDB group tree to discover SOMAs and populate the
     * SOMA URI map. This function is called recursively in order to discover
     * all SOMAs in all subgroups.
     *
     * @param group TileDB group to inspect for SOMAs and subgroups.
     * @param parent Hierarchical group name of the group's parent.
     */
    void build_uri_map(Group& group, std::string_view parent = "");
};

};  // namespace tiledbsoma

#endif
