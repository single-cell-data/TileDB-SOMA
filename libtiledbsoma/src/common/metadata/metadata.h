/**
 * @file   common.h
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
class MetadataCache {
   public:
    MetadataCache(tiledb::Array& array);
    MetadataCache(tiledb::Group& group);
    virtual ~MetadataCache() = default;

    bool contains(const std::string& key) const;

    std::map<std::string, MetadataValue> get() const;
    std::optional<MetadataValue> get(const std::string& key) const;

    void set(const std::string& key, MetadataValue value);

    void del(const std::string& key);

    void write(tiledb::Array& array);
    void write(tiledb::Group& group);

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