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

#include <tiledb/tiledb>

#include "enums.h"
#include "soma_dataframe.h"
#include "soma_dense_ndarray.h"
#include "soma_group.h"
#include "soma_object.h"
#include "soma_sparse_ndarray.h"

namespace tiledbsoma {

class SOMAExperiment;
class SOMAMeasurement;

using namespace tiledb;

class SOMACollection : public SOMAGroup {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMACollection object at the given URI.
     *
     * @param ctx TileDB context
     * @param uri URI to create the SOMACollection
     */
    static std::unique_ptr<SOMACollection> create(
        std::string_view uri, std::shared_ptr<SOMAContext> ctx);

    /**
     * @brief Open a group at the specified URI and return SOMACollection
     * object.
     *
     * @param uri URI of the array
     * @param mode read or write
     * @param ctx TileDB context
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::shared_ptr<SOMACollection> SOMACollection
     */
    static std::unique_ptr<SOMACollection> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMACollection object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param key key of the array
     * @param timestamp Optional pair indicating timestamp start and end
     */
    SOMACollection(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp)
        : SOMAGroup(
              mode,
              uri,
              ctx,
              std::string(std::filesystem::path(uri).filename()),  // group name
              timestamp){};

    SOMACollection(const SOMAGroup& other)
        : SOMAGroup(other) {
    }

    SOMACollection() = delete;
    SOMACollection(const SOMACollection&) = default;
    SOMACollection(SOMACollection&&) = default;
    ~SOMACollection() = default;

    using SOMAGroup::open;

    /**
     * Closes the SOMACollection object.
     */
    void close();

    /**
     * Get the SOMAObject associated with the key.
     *
     * @param key of member
     */
    std::shared_ptr<SOMAObject> get(const std::string& key);

    /**
     * Create and add a SOMACollection to the SOMACollection.
     *
     * @param key of collection
     * @param uri of SOMACollection to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     */
    std::shared_ptr<SOMACollection> add_new_collection(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx);

    /**
     * Create and add a SOMAExperiment to the SOMACollection.
     *
     * @param key of collection
     * @param uri of SOMAExperiment to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     */
    std::shared_ptr<SOMAExperiment> add_new_experiment(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        ArraySchema schema);

    /**
     * Create and add a SOMAMeasurement to the SOMACollection.
     *
     * @param key of collection
     * @param uri of SOMAMeasurement to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     */
    std::shared_ptr<SOMAMeasurement> add_new_measurement(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        ArraySchema schema);

    /**
     * Create and add a SOMADataFrame to the SOMACollection.
     *
     * @param key of dataframe
     * @param uri of SOMADataFrame to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     */
    std::shared_ptr<SOMADataFrame> add_new_dataframe(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        ArraySchema schema);

    /**
     * Create and add a SOMADenseNDArray to the SOMACollection.
     *
     * @param key of dense array
     * @param uri of SOMADenseNDArray to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     */
    std::shared_ptr<SOMADenseNDArray> add_new_dense_ndarray(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        ArraySchema schema);

    /**
     * Create and add a SOMASparseNDArray to the SOMACollection.
     *
     * @param key of sparse array
     * @param uri of SOMASparseNDArray to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     */
    std::shared_ptr<SOMASparseNDArray> add_new_sparse_ndarray(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<SOMAContext> ctx,
        ArraySchema schema);

   protected:
    //===================================================================
    //= protected non-static
    //===================================================================

    // Members of the SOMACollection
    std::map<std::string, std::shared_ptr<SOMAObject>> children_;
};
}  // namespace tiledbsoma

#endif  // SOMA_COLLECTION