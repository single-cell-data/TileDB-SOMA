/**
 * @file   common.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_COMMON_H
#define COMMON_COMMON_H

#include <stdexcept>
#include <string>
#include <string_view>
#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

namespace tiledbsoma::common {
/**
 * @brief Get a human redable name from the name of a type
 * 
 * @param name A mangled name as it comes from `std::type_info::name`
 * 
 * @return If compiled under GCC, returns the demangled type name, otherwise returns the mangled name as is.
 */
std::string demangle_name(std::string_view name);

/**
 * @brief Get the value associated with the given key from a TileDB config object.
 * 
 * @tparam T The datatype to cast the return value to
 * 
 * @param config The config object to extract the value from
 * @param key The key to fetch the requested value with
 * @param default_value The value to return if the requested key is missing
 * 
 * @return The value associated with the given key casted to the requested type `T` if present, otherwise return the default value
 */
template <typename T>
    requires std::floating_point<T> || std::integral<T> || std::same_as<T, std::string>
T get_config_value(const tiledb::Config& config, std::string_view key, T default_value) {
    if (!config.contains(key)) {
        return default_value;
    }

    std::string value_str = config.get(key.data());
    size_t pos;
    T value;

    if constexpr (std::is_floating_point_v<T>) {
        value = static_cast<T>(std::stold(value_str, &pos));
    } else if constexpr (std::is_integral_v<T>) {
        if constexpr (std::is_signed_v<T>) {
            value = static_cast<T>(std::stoll(value_str, &pos));
        } else {
            value = static_cast<T>(std::stoull(value_str, &pos));
        }
    } else {
        value = value_str;
        pos = value_str.size();
    }

    if (pos != value_str.size()) {
        throw std::runtime_error(
            "[Context][get_config_value] Unable to convert '" + std::string(key) + "' with value '" + value_str +
            "' to type '" + demangle_name(typeid(T).name()) + "'");
    }

    return value;
}

size_t enumeration_value_count(const tiledb::Context& ctx, const tiledb::Enumeration& enumeration);

/**
 * @brief Get the maximum number of Enumeration values the given index datatype can reference.
 * 
 * @param index_type The datatype of the index column.
 * @return The maximum number of enumeration values the index datatype can reference.
 */
size_t get_max_capacity(tiledb_datatype_t index_type);

};  // namespace tiledbsoma::common

#endif