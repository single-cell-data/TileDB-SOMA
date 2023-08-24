/**
 * @file   soma_group.h
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
 *   This declares the SOMAGroup class.
 */

#ifndef SOMA_GROUP
#define SOMA_GROUP

#include <future>
#include <stdexcept>
#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include "../utils/common.h"
#include "enums.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAGroup {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAGroup object at the given URI.
     *
     * @param ctx TileDB context
     * @param uri URI to create the SOMAGroup
     * @param soma_type SOMACollection, SOMAMeasurement, or SOMAExperiment
     */
    static void create(
        std::shared_ptr<Context> ctx,
        std::string_view uri,
        std::string soma_type);

    /**
     * @brief Open a group at the specified URI and return SOMAGroup
     * object.
     *
     * @param mode read or write
     * @param uri URI of the group
     * @param name Name of the group
     * @param platform_config Config parameter dictionary
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::unique_ptr<SOMAGroup> SOMAGroup
     */
    static std::unique_ptr<SOMAGroup> open(
        OpenMode mode,
        std::string_view uri,
        std::string_view name = "unnamed",
        std::map<std::string, std::string> platform_config = {},
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Open a group at the specified URI and return SOMAGroup
     * object.
     *
     * @param mode read or write
     * @param ctx TileDB context
     * @param uri URI of the group
     * @param name Name of the group
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::unique_ptr<SOMAGroup> SOMAGroup
     */
    static std::unique_ptr<SOMAGroup> open(
        OpenMode mode,
        std::shared_ptr<Context> ctx,
        std::string_view uri,
        std::string_view name = "unnamed",
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMAGroup object.
     *
     * @param mode read or write
     * @param uri URI of the group
     * @param name Name of the group
     * @param ctx TileDB context
     * @param timestamp Optional pair indicating timestamp start and end
     */
    SOMAGroup(
        OpenMode mode,
        std::string_view uri,
        std::string_view name,
        std::shared_ptr<Context> ctx,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    SOMAGroup() = delete;
    SOMAGroup(const SOMAGroup&) = delete;
    SOMAGroup(SOMAGroup&&) = default;
    ~SOMAGroup() = default;

    /**
     * Open the SOMAGroup object.
     *
     * @param mode read or write
     * @param timestamp Optional pair indicating timestamp start and end
     */
    void open(
        OpenMode mode,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * Close the SOMAGroup object.
     */
    void close();

    /**
     * Check if the SOMAGroup is open.
     *
     * @return bool true if open
     */
    bool is_open() const {
        return group_->is_open();
    }

    /**
     * Get the SOMAGroup URI.
     */
    const std::string uri() const;

    /**
     * Get the Context associated with the SOMAGroup.
     *
     * @return std::shared_ptr<Context>
     */
    std::shared_ptr<Context> ctx();

    /**
     * Get a member from the SOMAGroup given the index.
     *
     * @param index of member
     */
    tiledb::Object get_member(uint64_t index) const;

    /**
     * Get a member from the SOMAGroup given the name.
     *
     * @param name of member
     */
    tiledb::Object get_member(const std::string& name) const;

    /**
     * Check if the SOMAGroup contains the given name.
     *
     * @param name of member
     */
    bool has_member(const std::string& name);

    /**
     * Add a named member to a SOMAGroup.
     *
     * @param uri of member to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param name of member
     */
    void add_member(
        const std::string& uri, URIType uri_type, const std::string& name);

    /**
     * Get the number of members in the SOMAGroup.
     */
    uint64_t get_length() const;

    /**
     * Remove a named member from the SOMAGroup.
     *
     * @param name of member
     */
    void remove_member(const std::string& name);

    /**
     * Return a SOMAGroup member to URI mapping.
     *
     * @return std::map<std::string, std::string>
     */
    std::map<std::string, std::string> member_to_uri_mapping() const;

    /**
     * Set metadata key-value items to an open array. The array must
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
        const void* value);

    /**
     * Delete a metadata key-value item from an open group. The group must
     * be opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be deleted.
     *
     * @note The writes will take effect only upon closing the group.
     *
     * @note If the key does not exist, this will take no effect
     *     (i.e., the function will not error out).
     */
    void delete_metadata(const std::string& key);

    /**
     * @brief Given a key, get the associated value datatype, number of
     * values, and value in binary form. The group must be opened in READ mode,
     * otherwise the function will error out.
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
    std::map<std::string, MetadataValue> get_metadata();
    std::optional<MetadataValue> get_metadata(const std::string& key);

    /**
     * Check if the key exists in metadata from an open group. The group must
     * be opened in READ mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be checked. UTF-8 encodings
     *     are acceptable.
     * @return true if the key exists, else false.
     */
    bool has_metadata(const std::string& key);

    /**
     * Return then number of metadata items in an open group. The group must
     * be opened in READ mode, otherwise the function will error out.
     */
    uint64_t metadata_num() const;

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    /**
     * Fills the metadata and member-to-uri caches upon opening the array.
     */
    void fill_caches();

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMAGroup URI
    std::string uri_;

    // Name displayed in log messages
    std::string name_;

    // TileDBGroup associated with the SOMAGroup
    std::shared_ptr<Group> group_;

    // Metadata cache
    std::map<std::string, MetadataValue> metadata_;

    // Member-to-URI cache
    std::map<std::string, std::string> member_to_uri_;
};

}  // namespace tiledbsoma

#endif  // SOMA_GROUP
