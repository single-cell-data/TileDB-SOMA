/**
 * @file   util.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines utilities.
 */

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include <cstring>
#include "common/logging/impl/logger.h"
#include "utils/util.h"

namespace tiledbsoma::util {

template <typename T>
VarlenBufferPair to_varlen_buffers(std::vector<T> data, bool arrow) {
    size_t nbytes = 0;
    for (auto& elem : data) {
        nbytes += elem.size();
    }

    std::string result;
    std::vector<uint64_t> offsets(data.size() + 1);
    size_t offset = 0;
    size_t idx = 0;

    for (auto& elem : data) {
        result += elem;
        offsets[idx++] = offset;
        offset += elem.size();
    }
    offsets[idx] = offset;

    // Remove extra arrow offset when creating buffers for TileDB write
    if (!arrow) {
        offsets.pop_back();
    }

    return {result, offsets};
}

template VarlenBufferPair to_varlen_buffers(std::vector<std::string>, bool arrow);

MetadataEntry decode_metadata(common::DataType type, uint32_t elements, const void* data) {
    auto decode = [&]<typename T>() -> MetadataEntry {
        if constexpr (std::is_same_v<T, std::string>) {
            return MetadataEntry(std::string(static_cast<const char*>(data), elements));
        } else {
            if (elements == 1) {
                return MetadataEntry(reinterpret_cast<const T*>(data)[0]);
            } else {
                return MetadataEntry(
                    std::vector<T>(reinterpret_cast<const T*>(data), reinterpret_cast<const T*>(data) + elements));
            }
        }
    };

    switch (type) {
        case common::DataType::int8:
            return decode.template operator()<int8_t>();
        case common::DataType::int16:
            return decode.template operator()<int16_t>();
        case common::DataType::int32:
            return decode.template operator()<int32_t>();
        case common::DataType::int64:
            return decode.template operator()<int64_t>();
        case common::DataType::uint8:
            return decode.template operator()<uint8_t>();
        case common::DataType::uint16:
            return decode.template operator()<uint16_t>();
        case common::DataType::uint32:
            return decode.template operator()<uint32_t>();
        case common::DataType::uint64:
            return decode.template operator()<uint64_t>();
        case common::DataType::float32:
            return decode.template operator()<float>();
        case common::DataType::float64:
            return decode.template operator()<double>();
        case common::DataType::string_ascii:
        case common::DataType::string_utf8:
            return decode.template operator()<std::string>();
        default:
            throw std::runtime_error(fmt::format("Unsupported metadata type '{}'", common::getName(type)));
    }
}

bool is_tiledb_uri(std::string_view uri) {
    return uri.find("tiledb://") == 0;
}

std::string rstrip_uri(std::string_view uri) {
    return std::regex_replace(std::string(uri), std::regex("/+$"), "");
}

std::optional<std::vector<uint8_t>> bitmap_to_uint8(const uint8_t* bitmap, size_t length, size_t offset) {
    if (bitmap == nullptr) {
        return std::nullopt;
    }

    std::vector<uint8_t> casted(length);
    ArrowBitsUnpackInt8(bitmap, offset, length, reinterpret_cast<int8_t*>(casted.data()));
    return casted;
}

std::unique_ptr<uint8_t[]> bitmap_to_uint8_ptr(const uint8_t* bitmap, size_t length, size_t offset) {
    if (bitmap == nullptr) {
        return nullptr;
    }

    std::unique_ptr<uint8_t[]> bytemap = std::make_unique_for_overwrite<uint8_t[]>(length);
    ArrowBitsUnpackInt8(bitmap, offset, length, reinterpret_cast<int8_t*>(bytemap.get()));
    return bytemap;
}

std::shared_ptr<SOMAColumn> find_column_by_name(
    std::span<const std::shared_ptr<SOMAColumn>> columns, std::string_view name) {
    auto column_it = std::find_if(columns.begin(), columns.end(), [&](auto col) { return col->name() == name; });

    if (column_it == columns.end()) {
        throw TileDBSOMAError(fmt::format("Index column '{}' missing", name));
    }

    return *column_it;
}

std::string get_enmr_label(ArrowSchema* index_schema, ArrowSchema* value_schema) {
    std::string format(value_schema->format);
    format = (format == "u") ? "U" : (format == "z" ? "Z" : format);
    return std::string(index_schema->name) + "_" + format;
}

tiledb::Enumeration get_enumeration(
    std::shared_ptr<tiledb::Context> ctx,
    std::shared_ptr<tiledb::Array> arr,
    ArrowSchema* index_schema,
    ArrowSchema* value_schema) {
    std::string new_way = util::get_enmr_label(index_schema, value_schema);
    std::string old_way = std::string(index_schema->name);
    try {
        // New-style names of the form {attr_name}_{arrow_format}, e.g. "foo_U"
        // for attributes written by tiledbsoma >= 1.16.0
        return tiledb::ArrayExperimental::get_enumeration(*ctx, *arr, new_way);
    } catch (const std::exception& e) {
        // Old-style names of the form {attr_name}, e.g. "foo"
        // for attributes written by tiledbsoma < 1.16.0
        try {
            return tiledb::ArrayExperimental::get_enumeration(*ctx, *arr, old_way);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(
                fmt::format(
                    "[get_enumeration] Could not find enumeration with name '{} or "
                    "'{}'",
                    new_way,
                    old_way));
        }
    }
}

};  // namespace tiledbsoma::util
