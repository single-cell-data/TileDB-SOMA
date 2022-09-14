/**
 * @file   soma.h
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
 *   This declares the base soma object
 */

#ifndef SOMA_H
#define SOMA_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <mutex>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>
#include "tiledbsoma/soma_query.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAQuery;  // forward declaration

class SOMA {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open a SOMA at the specified URI and return a SOMA object.
     *
     * @param uri URI of SOMA
     * @param ctx TileDB context
     * @return SOMA object
     */
    static std::unique_ptr<SOMA> open(
        std::string_view uri,
        std::shared_ptr<Context> ctx = std::make_shared<Context>());

    /**
     * @brief Open a SOMA at the specified URI and return a SOMA object.
     *
     * @param uri URI of SOMA
     * @param config TileDB config
     * @return SOMA object
     */
    static std::unique_ptr<SOMA> open(
        std::string_view uri, const Config& config);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMA object
     *
     * @param uri URI of SOMA
     */
    SOMA(std::string_view uri, std::shared_ptr<Context> ctx);

    /**
     * @brief Return a map of hierarchical array names to array URIs for all
     * arrays in the SOMA.
     *
     * NOTE: If the SOMA URI is *not* a TileDB Cloud URI and the array URIs are
     * TileDB Cloud URIs, the array URIs are converted to relative URIs.
     *
     * @return std::unordered_map<std::string, std::string> Map of array name to
     * array URI
     */
    std::unordered_map<std::string, std::string> list_arrays();

    /**
     * @brief Open an array in the SOMA with the provided name, which is a
     * relative path. For example, for the "X array", name = "X/data".
     *
     * @param name Name of the array to open
     * @return Array TileDB array
     */
    std::shared_ptr<Array> open_array(const std::string& name);

    /**
     * @brief Create a SOMAQuery for this SOMA.
     *
     * @param name SOMA name
     * @return std::unique_ptr<SOMAQuery> A SOMA query
     */
    std::unique_ptr<SOMAQuery> query(std::string name = "") {
        return std::make_unique<SOMAQuery>(this, name);
    }

    /**
     * @brief Return TileDB context of SOMA.
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

    // SOMA URI
    std::string uri_;

    // Map of array name to array URI
    std::unordered_map<std::string, std::string> array_uri_map_;

    // Flag that is true if TileDB Cloud URIs were converted to relative URIs
    bool group_uri_override_ = false;

    // Mutex to control parallel access
    std::mutex mtx_;

    /**
     * @brief Walk the TileDB group tree to discover arrays and populate the
     * array URI map. This function is called recursively in order to discover
     * all arrays in all subgroups.
     *
     * @param group TileDB group to inspect for arrays and subgroups.
     * @param parent Hierarchical group name of the group's parent.
     */
    void build_uri_map(Group& group, std::string_view parent = "");
};

};  // namespace tiledbsoma

#endif
