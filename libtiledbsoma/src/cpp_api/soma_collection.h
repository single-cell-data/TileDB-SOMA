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
#include "soma_object.h"

namespace tiledbsoma {

class SOMAGroup;

using namespace tiledb;

class SOMACollection : SOMAObject {
   public:
    static std::unique_ptr<SOMACollection> open(
        std::shared_ptr<Context>, const std::string&);

    SOMACollection(std::shared_ptr<Context>, const std::string&);
    // std::map<std::string, std::string> platform_config);

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
     */
    SOMACollection get(const std::string& name);

    /**
     * Check if the SOMACollection contains the given name.
     */
    bool has(const std::string& name);

    /**
     * Get the number of SOMAObjects in the SOMACollection.
     */
    uint64_t length() const;

    /**
     * Delete the SOMAObject associated with the name.
     */
    void del(const std::string& name);

    /**
     * Get the member name to URI mapping of the SOMACollection.
     */
    std::map<std::string, std::string> member_to_uri_mapping() const;

    void add_new_collection(const std::string&, const SOMACollection&);
    // void add_new_dataframe(const std::string&, SOMADataFrame);
    // void add_new_dense_ndarray(const std::string&, SOMADenseNDArray);
    // void add_new_sparse_ndarray(const std::string&, SOMASparseNDArray);

   private:
    std::shared_ptr<Context> ctx_;
    std::unique_ptr<SOMAGroup> group_;
    std::string name_;
};
}  // namespace tiledbsoma

#endif  // SOMA_COLLECTION