/**
 * @file   metadata.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_METADATA_H
#define COMMON_METADATA_H

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "types.h"

#pragma region Forward declarations

namespace tiledb {
class Array;
class Group;
}  // namespace tiledb

#pragma endregion

namespace tiledbsoma::common {

/**
 * A metadata wrapper around a TileDB object providing a user friendly interface to access and act upon the metadata 
 * stored in the TileDB object.
 */
class MetadataCache {
   public:
    /**
     * @brief Construct a new metadata cache from TileDB Array metadata.
     * 
     * @param array The TileDB Array to read the metadata from. The array should be open in READ mode.
     */
    MetadataCache(tiledb::Array& array);

    /**
     * @brief Construct a new metadata cache from TileDB Group metadata.
     * 
     * @param array The TileDB Group to read the metadata from. The group should be open in READ mode.
     */
    MetadataCache(tiledb::Group& group);

    virtual ~MetadataCache() = default;

    /**
     * @brief Check if the cache contains the given key.
     * 
     * @param key The key to check if it exists.
     */
    bool contains(const std::string& key) const;

    /**
     * @brief Get a map containing all the existing typed metadata key value pairs.
     */
    std::map<std::string, MetadataValue> get() const;

    /**
     * @brief Get the typed value associated with the givem metadata key if it exists.
     * 
     * @param key The key to retrieve the value from.
     */
    std::optional<MetadataValue> get(const std::string& key) const;

    /**
     * @brief Set a typed metadata key value pair.
     * 
     * @param key The key to set the value to.
     * @param value The value to set.
     */
    void set(const std::string& key, MetadataValue value);

    /**
     * @brief Delete a metadata pair associated with the given key.
     * 
     * @param key The key to delete.
     */
    void del(const std::string& key);

    /**
     * @brief Write the changes of the cache back to the backing TileDB Array.
     * 
     * @param array The TileDB Array to write the changes to. The array should be open in WRITE mode.
     */
    void write(tiledb::Array& array);

    /**
     * @brief Write the changes of the cache back to the backing TileDB Group.
     * 
     * @param array The TileDB Group to write the changes to. The group should be open in WRITE mode.
     */
    void write(tiledb::Group& group);

    /**
     * @brief Get the number of keys in the cache.
     */
    std::size_t size() const;

   private:
    enum class DictMod { absent, added, present, updated, deleted };

    static DictMod next_state_(DictMod current_state, const std::string& action);

    DictMod current_state_(const std::string& key) const;

    std::map<std::string, MetadataValue> metadata_;
    std::map<std::string, DictMod> mods_;
};
}  // namespace tiledbsoma::common

#endif