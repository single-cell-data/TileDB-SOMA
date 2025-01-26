/**
 * @file   soma_group.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
#include "soma_object.h"

namespace tiledbsoma {
using namespace tiledb;

// Pair storing uri and soma type
using SOMAGroupEntry = std::pair<std::string, std::string>;

class SOMAGroup : public SOMAObject {
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
     * @param timestamp Optional pair indicating timestamp start and end
     */
    static std::unique_ptr<SOMAGroup> create(
        std::shared_ptr<SOMAContext> ctx,
        std::string_view uri,
        std::string_view soma_type,
        std::optional<TimestampRange> timestamp = std::nullopt);

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
        std::shared_ptr<SOMAContext> ctx,
        std::string_view name = "unnamed",
        std::optional<TimestampRange> timestamp = std::nullopt);

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
        std::shared_ptr<SOMAContext> ctx,
        std::string_view name,
        std::optional<TimestampRange> timestamp = std::nullopt);

    SOMAGroup(
        std::shared_ptr<SOMAContext> ctx,
        std::shared_ptr<Group> group,
        std::optional<TimestampRange> timestamp);

    SOMAGroup() = delete;
    SOMAGroup(const SOMAGroup&) = default;
    SOMAGroup(SOMAGroup&&) = default;
    virtual ~SOMAGroup() = default;

    /**
     * Open the SOMAGroup object.
     *
     * @param mode read or write
     * @param timestamp Optional pair indicating timestamp start and end
     */
    void open(
        OpenMode mode, std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Return a new SOMAGroup with the given mode at the current Unix timestamp.
     *
     * @param mode if the OpenMode is not given, If the SOMAObject was opened in
     * READ mode, reopen it in WRITE mode and vice versa
     * @param timestamp Timestamp
     */
    std::unique_ptr<SOMAGroup> reopen(
        OpenMode mode, std::optional<TimestampRange> timestamp = std::nullopt);

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
     * Get whether the SOMAGroup was open in read or write mode.
     *
     * @return OpenMode
     */
    OpenMode mode() const {
        return group_->query_type() == TILEDB_READ ? OpenMode::read :
                                                     OpenMode::write;
    }

    /**
     * Get the SOMAGroup URI.
     */
    const std::string uri() const;

    /**
     * Get the context associated with the SOMAGroup.
     *
     * @return SOMAContext
     */
    std::shared_ptr<SOMAContext> ctx();

    /**
     * Check if a named member is relative
     *
     * @param name of member to retrieve associated relative indicator.
     */
    bool is_relative(std::string name) const {
        return group_->is_relative(name);
    }

    /**
     * Get a member from the SOMAGroup given the index.
     *
     * @param index of member
     */
    tiledb::Object get(uint64_t index) const;

    /**
     * Get a member from the SOMAGroup given the name.
     *
     * @param name of member
     */
    tiledb::Object get(const std::string& name) const;

    /**
     * Check if the SOMAGroup contains the given name.
     *
     * @param name of member
     */
    bool has(const std::string& name);

    /**
     * Add a named member to a SOMAGroup.
     *
     * @param uri of member to add
     * @param uri_type whether the given URI is automatic (default), absolute,
     * or relative
     * @param name of member
     */
    void set(
        const std::string& uri,
        URIType uri_type,
        const std::string& name,
        const std::string& soma_type);

    /**
     * Get the number of members in the SOMAGroup.
     */
    uint64_t count() const;

    /**
     * Remove a named member from the SOMAGroup.
     *
     * @param name of member
     */
    void del(const std::string& name);

    /**
     * Return a mapping of all members in the group with its uri and type.
     *
     * @return std::optional<TimestampRange>
     */
    std::map<std::string, SOMAGroupEntry> members_map() const;

    /**
     * Return optional timestamp pair SOMAArray was opened with.
     */
    std::optional<TimestampRange> timestamp();

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
     * @param force A boolean toggle to suppress internal checks, defaults to
     *     false.
     *
     * @note The writes will take effect only upon closing the array.
     */
    void set_metadata(
        const std::string& key,
        tiledb_datatype_t value_type,
        uint32_t value_num,
        const void* value,
        bool force = false);

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
    void delete_metadata(const std::string& key, bool force = false);

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
     * Helper function to set the pass in timestamp in the config associated
     * with the SOMAContext passed in
     */
    static Config _set_timestamp(
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp);

    /**
     * Fills the metadata and member-to-uri caches upon opening the array.
     */
    void fill_caches();

    // SOMA context
    std::shared_ptr<SOMAContext> ctx_;

    // SOMAGroup URI
    std::string uri_;

    // Name displayed in log messages
    std::string name_;

    // TileDB Group associated with the SOMAGroup
    std::shared_ptr<Group> group_;

    // Metadata values need to be accessible in write mode as well. When adding
    // or deleting values in the group, instead of closing to update to
    // metadata; then reopening to read the group; and again reopening to
    // restore the group back to write mode, we just store the modifications to
    // this cache
    std::map<std::string, MetadataValue> metadata_;

    // Group associated with metadata_. We need to keep this read-mode group
    // alive in order for the metadata value pointers in the cache to be
    // accessible
    std::shared_ptr<Group> cache_group_;

    // Read timestamp range (start, end)
    std::optional<TimestampRange> timestamp_;

    // Member-to-URI cache
    std::map<std::string, SOMAGroupEntry> members_map_;
};

}  // namespace tiledbsoma

#endif  // SOMA_GROUP
