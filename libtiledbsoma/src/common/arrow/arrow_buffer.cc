/**
 * @file   arrow_buffer.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#include "arrow_buffer.h"
#include "nanoarrow/nanoarrow.hpp"

#include <format>
#include <limits>

namespace tiledbsoma::common::arrow {

size_t IArrowBufferStorage::length() const {
    return length_;
}

size_t IArrowBufferStorage::null_count() const {
    return null_count_;
}

ArrayArrowBufferStorage::ArrayArrowBufferStorage(tiledb_datatype_t type, size_t length, bool nullable) {
    length_ = length;
    if (type == TILEDB_BOOL) {
        data_size_ = (length + 7) / 8;
    } else {
        data_size_ = length * tiledb::impl::type_size(type);
    }

    data_buffer_ = std::make_unique_for_overwrite<std::byte[]>(data_size_);

    if (nullable) {
        validity_buffer_ = std::make_unique_for_overwrite<std::byte[]>((length + 7) / 8);
    } else {
        null_count_ = 0;
    }
}

ArrayArrowBufferStorage::ArrayArrowBufferStorage(
    tiledb_datatype_t offset_type, size_t length, size_t data_size, bool nullable) {
    length_ = length;
    data_size_ = data_size;
    offsets_size_ = (length + 1) * tiledb::impl::type_size(offset_type);

    data_buffer_ = std::make_unique_for_overwrite<std::byte[]>(data_size_);
    offset_buffer_ = std::make_unique_for_overwrite<std::byte[]>(offsets_size_);

    if (nullable) {
        validity_buffer_ = std::make_unique_for_overwrite<std::byte[]>((length + 7) / 8);
    } else {
        null_count_ = 0;
    }
}

ArrayArrowBufferStorage::ArrayArrowBufferStorage(
    tiledb_datatype_t type, size_t length, std::unique_ptr<std::byte[]> data, std::unique_ptr<uint8_t[]> validity) {
    length_ = length;

    if (type == TILEDB_BOOL) {
        data_size_ = (length + 7) / 8;
    } else {
        data_size_ = length * tiledb::impl::type_size(type);
    }

    if (validity) {
        null_count_ = length_ - ArrowBitCountSet(validity.get(), 0, length_);
    } else {
        null_count_ = 0;
    }

    data_buffer_ = std::move(data);
    validity_buffer_ = std::unique_ptr<std::byte[]>(reinterpret_cast<std::byte*>(validity.release()));
}

std::span<const std::byte> ArrayArrowBufferStorage::data() const {
    return std::span<const std::byte>(data_buffer_.get(), data_size_);
}

std::span<const std::byte> ArrayArrowBufferStorage::offsets() const {
    return std::span<const std::byte>(offset_buffer_.get(), offsets_size_);
}

std::span<const std::byte> ArrayArrowBufferStorage::validity() const {
    return std::span<const std::byte>(validity_buffer_.get(), (length_ + 7) / 8);
}

ArrowBuffer::ArrowBuffer(std::unique_ptr<IArrowBufferStorage> storage, std::string_view name)
    : storage_(std::move(storage))
    , name_(name) {
}

ArrowBuffer::ArrowBuffer(const tiledb::Enumeration& enumeration, bool large_offsets) {
    tiledb::Context ctx = enumeration.context();

    const void* offsets;
    uint64_t offsets_size;

    const void* data;
    uint64_t data_size;
    ctx.handle_error(tiledb_enumeration_get_data(ctx.ptr().get(), enumeration.ptr().get(), &data, &data_size));

    std::unique_ptr<std::byte[]> data_buffer = std::make_unique_for_overwrite<std::byte[]>(data_size);
    std::memcpy(data_buffer.get(), data, data_size);

    switch (enumeration.type()) {
        case TILEDB_CHAR:
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_BLOB:
        case TILEDB_GEOM_WKT:
        case TILEDB_GEOM_WKB: {
            ctx.handle_error(
                tiledb_enumeration_get_offsets(ctx.ptr().get(), enumeration.ptr().get(), &offsets, &offsets_size));

            size_t count = offsets_size / sizeof(uint64_t);

            if (large_offsets) {
                if (data_size > std::numeric_limits<int64_t>::max()) {
                    throw std::runtime_error(
                        std::format(
                            "[ArrowBuffer] Int64 cannot represent indices for `{}` enumeration: Datatype too small",
                            enumeration.name()));
                }

                std::unique_ptr<int64_t[]> offsets_buffer = std::make_unique_for_overwrite<int64_t[]>(count + 1);
                std::memcpy(offsets_buffer.get(), offsets, offsets_size);
                offsets_buffer[count] = static_cast<int64_t>(data_size);

                storage_ = std::make_unique<ArrayArrowBufferStorage>(
                    std::move(data_buffer), std::move(offsets_buffer), count);
            } else {
                if (data_size > std::numeric_limits<int32_t>::max()) {
                    throw std::runtime_error(
                        std::format(
                            "[ArrowBuffer] Int32 cannot represent indices for `{}` enumeration: Datatype too small",
                            enumeration.name()));
                }

                std::unique_ptr<int32_t[]> offsets_buffer = std::make_unique_for_overwrite<int32_t[]>(count + 1);
                std::span<const uint64_t> offsets_v(static_cast<const uint64_t*>(offsets), count);
                for (size_t i = 0; i < count; ++i) {
                    offsets_buffer[i] = static_cast<int32_t>(offsets_v[i]);
                }
                offsets_buffer[count] = static_cast<int32_t>(data_size);

                storage_ = std::make_unique<ArrayArrowBufferStorage>(
                    std::move(data_buffer), std::move(offsets_buffer), count);
            }
        } break;
        case TILEDB_BOOL: {
            std::span<const bool> data_v(static_cast<const bool*>(data), data_size);
            size_t count = data_size / sizeof(bool);

            // If the enumeration is not empty
            if (count > 0) {
                // Represent the Boolean vector with, at most, the last two
                // bits. In Arrow, Boolean values are LSB packed
                uint8_t packed_data = 0;
                for (size_t i = 0; i < count; ++i)
                    packed_data |= (data_v[i] << i);

                std::memcpy(data_buffer.get(), &packed_data, 1);
            }

            storage_ = std::make_unique<ArrayArrowBufferStorage>(enumeration.type(), count, std::move(data_buffer));
        } break;
        case TILEDB_INT8:
        case TILEDB_UINT8:
        case TILEDB_INT16:
        case TILEDB_UINT16:
        case TILEDB_INT32:
        case TILEDB_UINT32:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_INT64:
        case TILEDB_UINT64:
        case TILEDB_FLOAT32:
        case TILEDB_FLOAT64:
            storage_ = std::make_unique<ArrayArrowBufferStorage>(
                enumeration.type(), data_size / tiledb::impl::type_size(enumeration.type()), std::move(data_buffer));
            break;
        default:
            throw std::runtime_error(
                std::format(
                    "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                    tiledb::impl::type_to_str(enumeration.type())));
    }

    name_ = enumeration.name();
}

std::string ArrowBuffer::name() const {
    return name_;
}

IArrowBufferStorage* ArrowBuffer::storage() const {
    return storage_.get();
}
}  // namespace tiledbsoma::common::arrow