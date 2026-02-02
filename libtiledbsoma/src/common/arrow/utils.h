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

template <typename T>
std::tuple<std::unique_ptr<std::byte[]>, std::unique_ptr<uint64_t[]>, std::unique_ptr<uint8_t[]>> dictionary_to_values(
    ArrowSchema* schema, ArrowArray* array) {
    if (array->dictionary->n_buffers == 3) {
        throw std::runtime_error("[dictionary_to_values] Variable sized dictionaries are not supported");
    }

    std::span<const T> values(
        static_cast<const T*>(array->dictionary->buffers[1]) + array->dictionary->offset, array->dictionary->length);

    std::unique_ptr<std::byte[]> data_buffer = std::make_unique_for_overwrite<std::byte[]>(array->length * sizeof(T));
    std::span<T> data_view(reinterpret_cast<T*>(data_buffer.get()), array->length);

    auto extract_values = [&]<typename IndexType>() {
        std::span<const IndexType> indices(reinterpret_cast<const IndexType*>(array->buffers[1]), array->length);

        for (size_t i = 0; i < indices.size(); ++i) {
            data_view[i] = values[indices[i]];
        }
    };

    switch (to_tiledb_format(schema->format)) {
        case TILEDB_INT8:
            extract_values.template operator()<int8_t>();
            break;
        case TILEDB_UINT8:
            extract_values.template operator()<uint8_t>();
            break;
        case TILEDB_INT16:
            extract_values.template operator()<int16_t>();
            break;
        case TILEDB_UINT16:
            extract_values.template operator()<uint16_t>();
            break;
        case TILEDB_INT32:
            extract_values.template operator()<int32_t>();
            break;
        case TILEDB_UINT32:
            extract_values.template operator()<uint32_t>();
            break;
        case TILEDB_INT64:
            extract_values.template operator()<int64_t>();
            break;
        case TILEDB_UINT64:
            extract_values.template operator()<uint64_t>();
            break;
        default:
            throw std::runtime_error(
                "Saw invalid index type when trying to promote indexes to "
                "values");
    }

    std::unique_ptr<std::uint8_t[]> validity = nullptr;
    if (schema->flags & ARROW_FLAG_NULLABLE) {
        validity = bitmap_to_bytemap(static_cast<const uint8_t*>(array->buffers[0]), array->length, array->offset);
    }
}

template <>
std::tuple<std::unique_ptr<std::byte[]>, std::unique_ptr<uint64_t[]>, std::unique_ptr<uint8_t[]>>
dictionary_to_values<std::string>(ArrowSchema* schema, ArrowArray* array);

template <>
std::tuple<std::unique_ptr<std::byte[]>, std::unique_ptr<uint64_t[]>, std::unique_ptr<uint8_t[]>>
dictionary_to_values<bool>(ArrowSchema* schema, ArrowArray* array);

std::vector<std::string_view> dictionary_values_view(ArrowSchema* schema, ArrowArray* array);
}  // namespace tiledbsoma::common::arrow

#endif