/**
 * @file   soma_collection.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
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
 *   This file defines the SOMACollection class.
 */

#ifndef SOMA_COLLECTION
#define SOMA_COLLECTION

#include <map>
#include <memory>
#include <string>
#include <tiledb/tiledb>
#include "soma_dataframe.h"
#include "soma_dense_ndarray.h"
#include "soma_object.h"
#include "soma_sparse_ndarray.h"

namespace tiledbsoma {

class SOMAGroup;

using namespace tiledb;

class SOMACollection : SOMAObject {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open and return a SOMACollection object at the given URI.
     *
     * @param ctx TileDB context
     * @param name name of the array
     */
    static std::unique_ptr<SOMACollection> open(
        std::shared_ptr<Context> ctx, const std::string& uri);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMACollection object.
     *
     * @param ctx TileDB context
     * @param name name of the array
     */
    SOMACollection(std::shared_ptr<Context> ctx, const std::string& uri);
    // std::map<std::string, std::string> platform_config);

    /**
     * Returns the constant "SOMACollection".
     */
    std::string type() const {
        return "SOMACollection";
    }

    /**
     * Closes the SOMACollection object.
     */
    void close();

    /**
     * Get the SOMAGroup URI.
     */
    std::string uri() const;

    // void set(const std::string& name, SOMAObject& object, bool
    // use_relative_uri);

    /**
     * Get the SOMAObject associated with the name.
     *
     * @param name of member
     */
    SOMACollection get(const std::string& name);

    /**
     * Check if the SOMACollection contains the given name.
     *
     * @param name of member
     */
    bool has(const std::string& name);

    /**
     * Get the number of SOMAObjects in the SOMACollection.
     */
    uint64_t count() const;

    /**
     * Delete the SOMAObject associated with the name.
     *
     * @param name of member
     */
    void del(const std::string& name);

    /**
     * Get the member name to URI mapping of the SOMACollection.
     */
    std::map<std::string, std::string> member_to_uri_mapping() const;

    /**
     * Add a SOMACollection to the SOMACollection.
     *
     * @param name of collection
     * @param collection SOMACollection to add
     */
    void add_new_collection(
        const std::string& name, SOMACollection& collection);

    /**
     * Add a SOMADataFrame to the SOMACollection.
     *
     * @param name of dataframe
     * @param dataframe SOMADataFrame to add
     */
    void add_new_dataframe(const std::string& name, SOMADataFrame& dataframe);

    /**
     * Add a SOMADenseNDArray to the SOMACollection.
     *
     * @param name of dense array
     * @param dataframe SOMADenseNDArray to add
     */
    void add_new_dense_ndarray(
        const std::string& name, const SOMADenseNDArray& array);

    /**
     * Add a SOMASparseNDArray to the SOMACollection.
     *
     * @param name of sparse array
     * @param dataframe SOMASparseNDArray to add
     */
    void add_new_sparse_ndarray(
        const std::string& name, const SOMASparseNDArray& array);

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMAGroup
    std::unique_ptr<SOMAGroup> group_;
};
}  // namespace tiledbsoma

#endif  // SOMA_COLLECTION