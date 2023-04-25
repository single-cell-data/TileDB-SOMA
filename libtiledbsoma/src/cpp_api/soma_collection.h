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

#include "soma_object.h"

// #include "soma_dataframe.h"
// #include "soma_dense_ndarray.h"
// #include "soma_experiment.h"
// #include "soma_measurement.h"
// #include "soma_sparse_ndarray.h"
namespace tiledbsoma {
class SOMAGroup;
class SOMADataFrame;
class SOMAExperiment;
class SOMAMeasurement;
class SOMASparseNDArray;
class SOMADenseNDArray;

using namespace tiledb;

class SOMACollection : public SOMAObject {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open a group at the specified URI and return SOMACollection
     * object.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param uri URI of the array
     * @param platform_config Config parameter dictionary
     * @return std::shared_ptr<SOMACollection> SOMACollection
     */
    static std::shared_ptr<SOMACollection> open(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::map<std::string, std::string> platform_config = {});

    /**
     * @brief Open a group at the specified URI and return SOMACollection
     * object.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param ctx TileDB context
     * @param uri URI of the array
     * @return std::shared_ptr<SOMACollection> SOMACollection
     */
    static std::shared_ptr<SOMACollection> open(
        tiledb_query_type_t mode,
        std::shared_ptr<Context> ctx,
        std::string_view uri);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMACollection object.
     *
     * @param ctx TileDB context
     * @param key key of the array
     */
    SOMACollection(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::shared_ptr<Context> ctx);

    /**
     * Returns the constant "SOMACollection".
     */
    std::string type() const;

    /**
     * Closes the SOMACollection object.
     */
    void close();

    /**
     * Get the SOMACollection URI.
     */
    std::string uri() const;

    /**
     * Get the Context associated with the SOMACollection.
     *
     * @return std::shared_ptr<Context>
     */
    std::shared_ptr<Context> ctx();

    /**
     * Set an already existing SOMAObject with the given key.
     *
     * @param key of member
     * @param object SOMA object to add
     */
    void set(const std::string& key, SOMAObject& object);

    /**
     * Get the SOMAObject associated with the key.
     *
     * @param key of member
     */
    std::shared_ptr<SOMAObject> get(const std::string& key);

    /**
     * Check if the SOMACollection contains the given key.
     *
     * @param key of member
     */
    bool has(const std::string& key);

    /**
     * Get the number of SOMAObjects in the SOMACollection.
     */
    uint64_t count() const;

    /**
     * Delete the SOMAObject associated with the key.
     *
     * @param key of member
     */
    void del(const std::string& key);

    /**
     * Get the member key to URI mapping of the SOMACollection.
     */
    std::map<std::string, std::string> member_to_uri_mapping() const;

    /**
     * Create and add a SOMACollection to the SOMACollection.
     *
     * @param key of collection
     * @param uri of SOMACollection to add
     * @param relative whether the given URI is relative
     */
    SOMACollection add_new_collection(
        std::string_view key,
        std::string_view uri,
        bool relative,
        std::shared_ptr<Context> ctx);

    // /**
    //  * Create and add a SOMAExperiment to the SOMACollection.
    //  *
    //  * @param key of collection
    //  * @param uri of SOMAExperiment to add
    //  * @param relative whether the given URI is relative
    //  */
    // SOMAExperiment add_new_experiment(
    //     std::string_view key,
    //     std::string_view uri,
    //     bool relative,
    //     std::shared_ptr<Context> ctx,
    //     SOMADataFrame& obs,
    //     SOMACollection& ms);

    // /**
    //  * Create and add a SOMAMeasurement to the SOMACollection.
    //  *
    //  * @param key of collection
    //  * @param uri of SOMAMeasurement to add
    //  * @param relative whether the given URI is relative
    //  */
    // SOMAMeasurement add_new_measurement(
    //     std::string_view key,
    //     std::string_view uri,
    //     bool relative,
    //     std::shared_ptr<Context> ctx,
    //     SOMADataFrame& var,
    //     SOMACollection& X,
    //     SOMACollection& obsm,
    //     SOMACollection& obsp,
    //     SOMACollection& varm,
    //     SOMACollection& varp);

    /**
     * Create and add a SOMADataFrame to the SOMACollection.
     *
     * @param key of dataframe
     * @param uri of SOMADataFrame to add
     * @param relative whether the given URI is relative
     */
    SOMADataFrame add_new_dataframe(
        std::string_view key,
        std::string_view uri,
        bool relative,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        std::string_view result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp);

    /**
     * Create and add a SOMADenseNDArray to the SOMACollection.
     *
     * @param key of dense array
     * @param uri of SOMADenseNDArray to add
     * @param relative whether the given URI is relative
     */
    SOMADenseNDArray add_new_dense_ndarray(
        std::string_view key,
        std::string_view uri,
        bool relative,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        std::string_view result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp);

    /**
     * Create and add a SOMASparseNDArray to the SOMACollection.
     *
     * @param key of sparse array
     * @param uri of SOMASparseNDArray to add
     * @param relative whether the given URI is relative
     */
    SOMASparseNDArray add_new_sparse_ndarray(
        std::string_view key,
        std::string_view uri,
        bool relative,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        std::string_view result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp);

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMAGroup
    std::shared_ptr<SOMAGroup> group_;
};
}  // namespace tiledbsoma

#endif  // SOMA_COLLECTION