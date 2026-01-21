/**
 * @file   common.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This declares the column buffer API
 */

#ifndef COMMON_COMMON_H
#define COMMON_COMMON_H

#include <stdexcept>
#include <string>
#include <string_view>
#include <tiledb/tiledb>

namespace tiledbsoma::common {

enum class StatusCode { OK, ERROR };

class Status {
   public:
    Status();
    Status(StatusCode code, std::string_view origin, std::string_view message);
    ~Status() = default;

    StatusCode code() const {
        return code_;
    }

    std::string origin() const {
        return origin_;
    }

    std::string message() const {
        return message_;
    }

   private:
    StatusCode code_;
    std::string origin_;
    std::string message_;
};

template <typename T>
concept is_data_buffer = std::same_as<std::unique_ptr<std::byte[]>, T> ||
                         (std::is_pointer_v<T> &&
                          (std::same_as<std::remove_const_t<std::remove_pointer_t<T>>, void> ||
                           std::same_as<std::remove_const_t<std::remove_pointer_t<T>>, std::byte> ||
                           std::integral<std::remove_const_t<std::remove_pointer_t<T>>> ||
                           std::floating_point<std::remove_const_t<std::remove_pointer_t<T>>>));

template <typename T>
concept is_offset_buffer = std::same_as<T, std::unique_ptr<uint64_t[]>> ||
                           (std::is_pointer_v<T> &&
                            (std::same_as<std::remove_const_t<std::remove_pointer_t<T>>, uint32_t> ||
                             std::same_as<std::remove_const_t<std::remove_pointer_t<T>>, uint64_t>));

std::string demangle_name(std::string_view name);

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

size_t get_max_capacity(tiledb_datatype_t index_type);

};  // namespace tiledbsoma::common

#endif