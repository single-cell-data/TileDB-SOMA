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

class SOMACollection : public SOMAObject {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMACollection object at the given URI.
     *
     * @param uri URI of the array
     * @param platform_config Optional config parameter dictionary
     */
    static std::unique_ptr<SOMACollection> create(
        std::string_view uri,
        std::map<std::string, std::string> platform_config = {});

    /**
     * @brief Create a SOMACollection object at the given URI.
     *
     * @param ctx TileDB context
     * @param uri URI to create the SOMACollection
     */
    static std::unique_ptr<SOMACollection> create(
        std::string_view uri, std::shared_ptr<Context> ctx);

    /**
     * @brief Open a group at the specified URI and return SOMACollection
     * object.
     *
     * @param uri URI of the array
     * @param mode read or write
     * @param platform_config Config parameter dictionary
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::shared_ptr<SOMACollection> SOMACollection
     */
    static std::unique_ptr<SOMACollection> open(
        std::string_view uri,
        OpenMode mode,
        std::map<std::string, std::string> platform_config = {},
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

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
        std::shared_ptr<Context> ctx,
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
        std::shared_ptr<Context> ctx,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp);

    SOMACollection() = delete;
    SOMACollection(const SOMACollection&) = default;
    SOMACollection(SOMACollection&&) = default;
    ~SOMACollection() = default;

    /**
     * Open the SOMACollection object.
     *
     * @param mode read or write
     * @param timestamp Timestamp
     */
    void open(
        OpenMode mode, std::optional<std::pair<uint64_t, uint64_t>> timestamp);

    /**
     * Closes the SOMACollection object.
     */
    void close();

    /**
     * Check if the SOMACollection is open.
     *
     * @return bool true if open
     */
    bool is_open() const {
        return group_->is_open();
    }

    /**
     * Return the constant "SOMACollection".
     *
     * @return std::string
     */
    const std::string type() const {
        return "SOMACollection";
    }

    /**
     * Get the SOMACollection URI.
     */
    const std::string uri() const;

    /**
     * Get the Context associated with the SOMACollection.
     *
     * @return std::shared_ptr<Context>
     */
    std::shared_ptr<Context> ctx();

    /**
     * Set an already existing SOMAObject at uri to the given key.
     *
     * @param uri of member to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param key to add
     */
    void set(std::string_view uri, URIType uri_type, const std::string& key);

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
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     */
    std::shared_ptr<SOMACollection> add_new_collection(
        std::string_view key,
        std::string_view uri,
        URIType uri_type,
        std::shared_ptr<Context> ctx);

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
        std::shared_ptr<Context> ctx,
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
        std::shared_ptr<Context> ctx,
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
        std::shared_ptr<Context> ctx,
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
        std::shared_ptr<Context> ctx,
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
        std::shared_ptr<Context> ctx,
        ArraySchema schema);

    /**
     * Set metadata key-value items to a SOMACollection. The SOMACollection must
     * opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be added. UTF-8 encodings
     *     are acceptable.
     * @param value_type The datatype of the value.
     * @param value_num The value may consist of more than one items of the
     *     same datatype. This argument indicates the number of items in the
     *     value component of the metadata.
     * @param value The metadata value in binary form.
     *
     * @note The writes will take effect only upon closing the array.
     */
    void set_metadata(
        const std::string& key,
        tiledb_datatype_t value_type,
        uint32_t value_num,
        const void* value) {
        group_->set_metadata(key, value_type, value_num, value);
    }

    /**
     * Delete a metadata key-value item from an open SOMACollection. The
     * SOMACollection must be opened in WRITE mode, otherwise the function will
     * error out.
     *
     * @param key The key of the metadata item to be deleted.
     *
     * @note The writes will take effect only upon closing the group.
     *
     * @note If the key does not exist, this will take no effect
     *     (i.e., the function will not error out).
     */
    void delete_metadata(const std::string& key) {
        group_->delete_metadata(key);
    }

    /**
     * @brief Given a key, get the associated value datatype, number of
     * values, and value in binary form.
     *
     * The value may consist of more than one items of the same datatype. Keys
     * that do not exist in the metadata will be return NULL for the value.
     *
     * **Example:**
     * @code{.cpp}
     * // Open the group for reading
     * tiledbsoma::SOMAGroup soma_group = SOMAGroup::open(TILEDB_READ,
     "s3://bucket-name/group-name");
     * tiledbsoma::MetadataValue meta_val = soma_group->get_metadata("key");
     * std::string key = std::get<MetadataInfo::key>(meta_val);
     * tiledb_datatype_t dtype = std::get<MetadataInfo::dtype>(meta_val);
     * uint32_t num = std::get<MetadataInfo::num>(meta_val);
     * const void* value = *((const
     int32_t*)std::get<MetadataInfo::value>(meta_val));
     * @endcode
     *
     * @param key The key of the metadata item to be retrieved. UTF-8 encodings
     *     are acceptable.
     * @return MetadataValue (std::tuple<std::string, tiledb_datatype_t,
     * uint32_t, const void*>)
     */
    std::optional<MetadataValue> get_metadata(const std::string& key) {
        return group_->get_metadata(key);
    }

    /**
     * Get a mapping of all metadata keys with its associated value datatype,
     * number of values, and value in binary form.
     *
     * @return std::map<std::string, MetadataValue>
     */
    std::map<std::string, MetadataValue> get_metadata() {
        return group_->get_metadata();
    }

    /**
     * Check if the key exists in metadata from an open SOMACollection.
     *
     * @param key The key of the metadata item to be checked. UTF-8 encodings
     *     are acceptable.
     * @return true if the key exists, else false.
     */
    bool has_metadata(const std::string& key) {
        return group_->has_metadata(key);
    }

    /**
     * Return then number of metadata items in an open SOMACollection. The group
     * must be opened in READ mode, otherwise the function will error out.
     */
    uint64_t metadata_num() const {
        return group_->metadata_num();
    }

   protected:
    //===================================================================
    //= protected non-static
    //===================================================================

    // SOMAGroup
    std::shared_ptr<SOMAGroup> group_;

    // Members of the SOMACollection
    std::map<std::string, std::shared_ptr<SOMAObject>> children_;
};
}  // namespace tiledbsoma

#endif  // SOMA_COLLECTION