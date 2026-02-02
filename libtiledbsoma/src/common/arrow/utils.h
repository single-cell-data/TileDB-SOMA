/**
 * @file   utils.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_ARROW_UTILS_H
#define COMMON_ARROW_UTILS_H

#include <tiledb/tiledb>
#include <unordered_map>

#include <span>
#include <string_view>
#include <tuple>

#include "nanoarrow/nanoarrow.hpp"

namespace tiledbsoma::common::arrow {
template <typename T>
using managed_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;
using ArrowTable = std::pair<managed_unique_ptr<ArrowArray>, managed_unique_ptr<ArrowSchema>>;

template <typename T, typename... Args>
    requires std::same_as<T, ArrowArray> || std::same_as<T, ArrowSchema>
managed_unique_ptr<T> make_managed_unique(Args&&... args) {
    return managed_unique_ptr<T>(new T(std::forward<Args>(args)...), [](T* arrow_struct) {
        if (arrow_struct->release != nullptr) {
            arrow_struct->release(arrow_struct);
        }

        delete arrow_struct;
    });
}

ArrowTable make_empty_arrow_table(std::string_view name, std::string_view format, size_t num_children);

tiledb_datatype_t to_tiledb_format(std::string_view format, std::string_view dtype_metadata = {});
std::string_view to_arrow_format(tiledb_datatype_t format, bool use_large = true);
std::string_view to_arrow_readable(std::string_view format);

std::unique_ptr<uint8_t[]> bitmap_to_bytemap(const uint8_t* bitmap, size_t length, size_t offset);
std::unique_ptr<uint8_t[]> bytemap_to_bitmap(const uint8_t* bytemap, size_t length, size_t offset);
}  // namespace tiledbsoma::common::arrow

#endif