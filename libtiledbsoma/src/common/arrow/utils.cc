/**
 * @file   utils.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#include <numeric>
#include <stdexcept>
#include <unordered_map>

#include "../logging/impl/logger.h"
#include "arrow_buffer.h"
#include "utils.h"

namespace tiledbsoma::common::arrow {

namespace impl {
constexpr std::array<char, 2> bool_dict({{0, 1}});
}

ArrowTable make_empty_arrow_table(std::string_view name, std::string_view format, size_t num_children = 0) {
    managed_unique_ptr<ArrowSchema> schema = make_managed_unique<ArrowSchema>();
    managed_unique_ptr<ArrowArray> array = make_managed_unique<ArrowArray>();

    // Copy the name string view contents to a null terminated char buffer
    schema->name = new char[name.size() + 1]();
    std::memcpy(const_cast<char*>(schema->name), name.data(), name.size());

    // Copy the format string view contents to a null terminated char buffer
    schema->format = new char[format.size() + 1]();
    std::memcpy(const_cast<char*>(schema->format), format.data(), format.size());

    schema->children = new ArrowSchema*[num_children];
    schema->dictionary = nullptr;
    schema->flags = 0;
    schema->metadata = nullptr;
    schema->n_children = static_cast<int64_t>(num_children);
    schema->private_data = nullptr;
    schema->release = [](ArrowSchema* schema) {
        delete[] schema->name;
        delete[] schema->format;
        delete[] schema->metadata;

        if (schema->children) {
            for (int64_t i = 0; i < schema->n_children; ++i) {
                schema->children[i]->release(schema->children[i]);
                delete schema->children[i];
            }

            delete[] schema->children;
        }

        if (schema->dictionary) {
            schema->dictionary->release(schema->dictionary);
            delete schema->dictionary;
        }

        schema->release = nullptr;
    };

    array->buffers = const_cast<const void**>(new void*[3]());
    array->children = new ArrowArray*[num_children];
    array->dictionary = nullptr;
    array->length = 0;
    array->n_buffers = 0;
    array->n_children = num_children;
    array->null_count = 0;
    array->offset = 0;
    array->private_data = new PrivateArrowBuffer();
    array->release = [](ArrowArray* array) {
        delete[] array->buffers;

        delete reinterpret_cast<PrivateArrowBuffer*>(array->private_data);

        if (array->children) {
            for (int64_t i = 0; i < array->n_children; ++i) {
                array->children[i]->release(array->children[i]);
                delete array->children[i];
            }

            delete[] array->children;
        }

        if (array->dictionary) {
            array->dictionary->release(array->dictionary);
            delete array->dictionary;
        }

        array->release = nullptr;
    };

    return std::make_pair(std::move(array), std::move(schema));
}

tiledb_datatype_t to_tiledb_format(std::string_view format, std::string_view dtype_metadata) {
    std::map<std::string_view, tiledb_datatype_t> _to_tiledb_format_map = {
        {"u", TILEDB_STRING_UTF8},    {"U", TILEDB_STRING_UTF8},
        {"z", TILEDB_CHAR},           {"Z", TILEDB_CHAR},
        {"c", TILEDB_INT8},           {"C", TILEDB_UINT8},
        {"s", TILEDB_INT16},          {"S", TILEDB_UINT16},
        {"i", TILEDB_INT32},          {"I", TILEDB_UINT32},
        {"l", TILEDB_INT64},          {"L", TILEDB_UINT64},
        {"f", TILEDB_FLOAT32},        {"g", TILEDB_FLOAT64},
        {"b", TILEDB_BOOL},           {"tss:", TILEDB_DATETIME_SEC},
        {"tsm:", TILEDB_DATETIME_MS}, {"tsu:", TILEDB_DATETIME_US},
        {"tsn:", TILEDB_DATETIME_NS},
    };

    try {
        auto dtype = _to_tiledb_format_map.at(format);

        if (dtype == TILEDB_BLOB && dtype_metadata.compare("WKB") == 0) {
            dtype = TILEDB_GEOM_WKB;
        } else if (dtype == TILEDB_STRING_UTF8 && dtype_metadata.compare("WKT") == 0) {
            dtype = TILEDB_GEOM_WKT;
        }

        return dtype;
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(
            fmt::format("ArrowAdapter: Unsupported Arrow type: {} ({})", format, to_arrow_readable(format)));
    }
}

std::string_view to_arrow_format(tiledb_datatype_t format, bool use_large) {
    auto u = use_large ? "U" : "u";
    auto z = use_large ? "Z" : "z";
    std::map<tiledb_datatype_t, std::string_view> _to_arrow_format_map = {
        {TILEDB_STRING_ASCII, u},     {TILEDB_CHAR, z},
        {TILEDB_STRING_UTF8, u},      {TILEDB_BLOB, z},
        {TILEDB_INT8, "c"},           {TILEDB_UINT8, "C"},
        {TILEDB_INT16, "s"},          {TILEDB_UINT16, "S"},
        {TILEDB_INT32, "i"},          {TILEDB_UINT32, "I"},
        {TILEDB_INT64, "l"},          {TILEDB_UINT64, "L"},
        {TILEDB_FLOAT32, "f"},        {TILEDB_FLOAT64, "g"},
        {TILEDB_BOOL, "b"},           {TILEDB_DATETIME_SEC, "tss:"},
        {TILEDB_DATETIME_MS, "tsm:"}, {TILEDB_DATETIME_US, "tsu:"},
        {TILEDB_DATETIME_NS, "tsn:"}, {TILEDB_GEOM_WKB, z},
        {TILEDB_GEOM_WKT, u}};
    try {
        return _to_arrow_format_map.at(format);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(
            fmt::format("ArrowAdapter: Unsupported TileDB type: {} ", tiledb::impl::type_to_str(format)));
    }
}

std::string_view to_arrow_readable(std::string_view format) {
    std::map<std::string_view, std::string_view> _to_arrow_readable = {
        {"n", "null"},
        {"b", "boolean"},
        {"c", "int8"},
        {"C", "uint8"},
        {"s", "int16"},
        {"S", "uint16"},
        {"i", "int32"},
        {"I", "uint32"},
        {"l", "int64"},
        {"L", "uint64"},
        {"e", "float16"},
        {"f", "float32"},
        {"g", "float64"},
        {"z", "binary"},
        {"Z", "large binary"},
        {"vz", "binary view"},
        {"u", "utf-8 string"},
        {"U", "large utf-8 string"},
        {"vu", "utf-8 view"},
        {"tdD", "date32 [days]"},
        {"tdm", "date64 [milliseconds]"},
        {"tts", "time32 [seconds]"},
        {"ttm", "time32 [milliseconds]"},
        {"ttu", "time64 [microseconds]"},
        {"ttn", "time64 [nanoseconds]"},
        {"tDs", "duration [seconds]"},
        {"tDm", "duration [milliseconds]"},
        {"tDu", "duration [microseconds]"},
        {"tDn", "duration [nanoseconds]"},
        {"tiM", "interval [months]"},
        {"tiD", "interval [days, time]"},
        {"tin", "interval [month, day, nanoseconds]"},
        {"+l", "list"},
        {"+L", "large list"},
        {"+vl", "list-view"},
        {"+vL", "large list-view"},
        {"+s", "struct"},
        {"+m", "map"},
        {"+r", "run-end encoded"}};

    auto it = _to_arrow_readable.find(format);
    return it != _to_arrow_readable.end() ? it->second :
                                            "unknown Arrow type [see "
                                            "https://arrow.apache.org/docs/format/"
                                            "CDataInterface.html#data-type-description-format-strings]";
}

std::unique_ptr<uint8_t[]> bitmap_to_bytemap(const uint8_t* bitmap, size_t length, size_t offset) {
    if (bitmap == nullptr) {
        return nullptr;
    }

    std::unique_ptr<uint8_t[]> bytemap = std::make_unique_for_overwrite<uint8_t[]>(length);
    ArrowBitsUnpackInt8(bitmap, offset, length, reinterpret_cast<int8_t*>(bytemap.get()));
    return bytemap;
}

std::unique_ptr<uint8_t[]> bytemap_to_bitmap(const uint8_t* bytemap, size_t length, size_t offset) {
    if (bytemap == nullptr) {
        return nullptr;
    }

    std::unique_ptr<uint8_t[]> bitmap = std::make_unique_for_overwrite<uint8_t[]>((length + 7) / 8);
    size_t i_dst = 0;
    for (size_t i_src = 0; i_src < length - offset; ++i_src) {
        // Overwrite every 8 bytes with a one-byte bitmap
        if (i_src % 8 == 0) {
            // Each bit in the bitmap corresponds to one byte in the bytemap
            // Note: the bitmap must be byte-aligned (8 bits)
            bitmap[i_dst] = 0;
            for (size_t i = i_src; i < i_src + 8 && i < length - offset; i++) {
                bitmap[i_dst] |= bytemap[i + offset] << (i % 8);
            }
            ++i_dst;
        }
    }

    return bitmap;
}

template <>
std::tuple<std::unique_ptr<std::byte[]>, std::unique_ptr<uint64_t[]>, std::unique_ptr<uint8_t[]>>
dictionary_to_values<std::string>(ArrowSchema* schema, ArrowArray* array) {
    auto extract_values = []<typename OffsetType>(ArrowSchema* schema, ArrowArray* array) {
        auto value_array = array->dictionary;

        std::span<const OffsetType> value_offsets(
            (OffsetType*)value_array->buffers[1] + value_array->offset, value_array->length);
        std::span<const std::byte> value_data((std::byte*)value_array->buffers[2], value_offsets.back());

        std::unique_ptr<uint64_t[]> offsets = std::make_unique_for_overwrite<uint64_t[]>(array->length + 1);
        std::unique_ptr<std::byte[]> data;

        auto calculate_value_size = [&]<typename S>() {
            // Calculate the total number of bytes to allocate for the values by summing the differences of the offset
            std::span<const S> indices(reinterpret_cast<const S*>(array->buffers[1]) + array->offset, array->length);

            size_t data_buffer_size = std::transform_reduce(
                indices.begin(), indices.end(), 0L, std::plus{}, [&](const auto& index) {
                    return value_offsets[index + 1] - value_offsets[index];
                });

            data = std::make_unique_for_overwrite<std::byte[]>(data_buffer_size);

            size_t offset = 0;
            offsets[0] = 0;

            for (size_t i = 0; i < indices.size(); ++i) {
                auto index = indices[i];
                std::memcpy(
                    data.get() + offset,
                    value_data.data() + value_offsets[index],
                    value_offsets[index + 1] - value_offsets[index]);

                offset += value_offsets[index + 1] - value_offsets[index];
                offsets[i + 1] = offset;
            }
        };

        switch (to_tiledb_format(schema->format)) {
            case TILEDB_INT8:
                calculate_value_size.template operator()<int8_t>();
                break;
            case TILEDB_UINT8:
                calculate_value_size.template operator()<uint8_t>();
                break;
            case TILEDB_INT16:
                calculate_value_size.template operator()<int16_t>();
                break;
            case TILEDB_UINT16:
                calculate_value_size.template operator()<uint16_t>();
                break;
            case TILEDB_INT32:
                calculate_value_size.template operator()<int32_t>();
                break;
            case TILEDB_UINT32:
                calculate_value_size.template operator()<uint32_t>();
                break;
            case TILEDB_INT64:
                calculate_value_size.template operator()<int64_t>();
                break;
            case TILEDB_UINT64:
                calculate_value_size.template operator()<uint64_t>();
                break;
            default:
                throw std::runtime_error("Saw invalid index type when trying to promote indexes to values");
        }

        return std::make_pair<std::unique_ptr<std::byte[]>, std::unique_ptr<uint64_t[]>>(
            std::move(data), std::move(offsets));
    };

    std::unique_ptr<std::byte[]> data_buffer;
    std::unique_ptr<uint64_t[]> offset_buffer;

    if ((strcmp(schema->dictionary->format, "U") == 0) || (strcmp(schema->dictionary->format, "Z") == 0)) {
        std::tie(data_buffer, offset_buffer) = extract_values.operator()<uint64_t>(schema, array);
    } else {
        std::tie(data_buffer, offset_buffer) = extract_values.operator()<uint32_t>(schema, array);
    }

    std::unique_ptr<std::uint8_t[]> validity = nullptr;
    if (schema->flags & ARROW_FLAG_NULLABLE) {
        validity = bitmap_to_bytemap(static_cast<const uint8_t*>(array->buffers[0]), array->length, array->offset);
    }

    return std::make_tuple(std::move(data_buffer), std::move(offset_buffer), std::move(validity));
}

template <>
std::tuple<std::unique_ptr<std::byte[]>, std::unique_ptr<uint64_t[]>, std::unique_ptr<uint8_t[]>>
dictionary_to_values<bool>(ArrowSchema* schema, ArrowArray* array) {
    std::unique_ptr<std::byte[]> data_buffer = std::make_unique_for_overwrite<std::byte[]>(array->length);
    std::span<uint8_t> data_view(reinterpret_cast<uint8_t*>(data_buffer.get()), array->length);

    auto extract_values = [&]<typename T>() {
        std::span<const T> indices(reinterpret_cast<const T*>(array->buffers[1]) + array->offset, array->length);

        for (int64_t i = 0; i < array->length; ++i) {
            data_view[i] = static_cast<uint8_t>(
                ArrowBitGet((const uint8_t*)array->dictionary->buffers[1], indices[i] + array->dictionary->offset));
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
            throw std::runtime_error("Saw invalid index type when trying to promote indexes to values");
    }

    std::unique_ptr<std::uint8_t[]> validity = nullptr;
    if (schema->flags & ARROW_FLAG_NULLABLE) {
        validity = bitmap_to_bytemap(static_cast<const uint8_t*>(array->buffers[0]), array->length, array->offset);
    }

    return std::make_tuple(std::move(data_buffer), std::unique_ptr<uint64_t[]>(nullptr), std::move(validity));
}

std::vector<std::string_view> dictionary_values_view(ArrowSchema* schema, ArrowArray* array) {
    std::vector<std::string_view> result;

    auto extract_var_size_view = [&]<typename OffsetType>() {
        if (array->length == 0) {
            return;
        }

        std::span<const OffsetType> offsets(
            static_cast<const OffsetType*>(array->buffers[1]) + array->offset, array->length + 1);
        std::string_view data(static_cast<const char*>(array->buffers[2]), offsets.back());

        result.reserve(array->length);

        for (size_t i = 0; i < offsets.size() - 1; i++) {
            const size_t length = offsets[i + 1] - offsets[i];
            result.push_back(data.substr(offsets[i], length));
        }
    };

    if (array->n_buffers == 3) {
        if (strcmp(schema->format, "u") == 0 || strcmp(schema->format, "z") == 0) {
            extract_var_size_view.template operator()<int32_t>();
        } else if (strcmp(schema->format, "U") == 0 || strcmp(schema->format, "Z") == 0) {
            extract_var_size_view.template operator()<int64_t>();
        } else {
            throw std::runtime_error(
                fmt::format("[dictionary_values_view] Unknown format '{}' for var sized dictionary.", schema->format));
        }
    } else {
        result.reserve(array->length);

        if (to_tiledb_format(schema->format) == TILEDB_BOOL) {
            std::unique_ptr<uint8_t[]> byteset = bitmap_to_bytemap(
                static_cast<const uint8_t*>(array->buffers[1]), array->length, array->offset);

            for (size_t i = 0; i < static_cast<size_t>(array->length); ++i) {
                if (byteset[i] == 0) {
                    result.emplace_back(impl::bool_dict.begin(), 1);
                } else {
                    result.emplace_back(impl::bool_dict.begin() + 1, 1);
                }
            }
        } else {
            const size_t stride = tiledb::impl::type_size(to_tiledb_format(schema->format));
            std::string_view data(
                static_cast<const char*>(array->buffers[1]) + array->offset * stride, array->length * stride);

            for (int64_t i = 0; i < array->length; ++i) {
                result.push_back(data.substr(i * stride, stride));
            }
        }
    }

    return result;
}
}  // namespace tiledbsoma::common::arrow
