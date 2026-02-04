/**
 * @file   arrow_buffer.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_ARROW_BUFFER_H
#define COMMON_ARROW_BUFFER_H

#include <concepts>
#include <limits>
#include <span>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include "../common.h"
#include "interface.h"
#include "nanoarrow/nanoarrow.hpp"

namespace tiledbsoma::common::arrow {

class ArrayArrowBufferStorage : public IArrowBufferStorage {
   public:
    ArrayArrowBufferStorage(tiledb_datatype_t type, size_t length, bool nullable = false);
    ArrayArrowBufferStorage(tiledb_datatype_t offset_type, size_t length, size_t data_size, bool nullable = false);
    ArrayArrowBufferStorage(
        tiledb_datatype_t type,
        size_t length,
        std::unique_ptr<std::byte[]> data,
        std::unique_ptr<uint8_t[]> validity = nullptr);

    template <typename IndexType>
        requires std::same_as<IndexType, uint64_t> || std::same_as<IndexType, int64_t> ||
                 std::same_as<IndexType, int32_t>
    ArrayArrowBufferStorage(
        std::unique_ptr<std::byte[]> data,
        std::unique_ptr<IndexType[]> offsets,
        size_t length,
        std::unique_ptr<uint8_t[]> validity = nullptr) {
        length_ = length;
        data_size_ = static_cast<size_t>(offsets[length_]);
        offsets_size_ = (length_ + 1) * sizeof(IndexType);

        if (data_size_ > std::numeric_limits<int64_t>::max()) {
            throw std::runtime_error(
                "Variable size buffer contains values unable to be represented by 64bit signed integers");
        }

        if (validity) {
            null_count_ = length_ - ArrowBitCountSet(validity.get(), 0, length_);
        } else {
            null_count_ = 0;
        }

        data_buffer_ = std::move(data);
        offset_buffer_ = std::unique_ptr<std::byte[]>(reinterpret_cast<std::byte*>(offsets.release()));
        validity_buffer_ = std::unique_ptr<std::byte[]>(reinterpret_cast<std::byte*>(validity.release()));
    }

    std::span<const std::byte> data() const override;
    std::span<const std::byte> offsets() const override;
    std::span<const std::byte> validity() const override;

   private:
    std::unique_ptr<std::byte[]> data_buffer_;
    std::unique_ptr<std::byte[]> offset_buffer_;
    std::unique_ptr<std::byte[]> validity_buffer_;
};

class ArrowBuffer {
   public:
    ArrowBuffer(std::unique_ptr<IArrowBufferStorage> storage, std::string_view name);
    ArrowBuffer(const tiledb::Enumeration& enumeration, bool large_offsets = true);

    IArrowBufferStorage* storage() const;
    std::string name() const;

   private:
    std::unique_ptr<IArrowBufferStorage> storage_;
    std::string name_;
};

struct PrivateArrowBuffer {
    PrivateArrowBuffer(const std::shared_ptr<ArrowBuffer>& buffer = nullptr)
        : buffer_(buffer) {
    }

    std::shared_ptr<ArrowBuffer> buffer_;
};
}  // namespace tiledbsoma::common::arrow

#endif