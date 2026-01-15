/**
 * @file   column_buffer.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the a ColumBuffer class.
 */

#include "column_buffer.h"
#include "common/logging/impl/logger.h"
#include "common/logging/logger.h"

namespace tiledbsoma {

using namespace tiledb;
using namespace common::logging;

#pragma region ColumnBuffer

#pragma region public static

void ColumnBuffer::to_bitmap(std::span<uint8_t> bytemap) {
    int i_dst = 0;
    for (unsigned int i_src = 0; i_src < bytemap.size(); i_src++) {
        // Overwrite every 8 bytes with a one-byte bitmap
        if (i_src % 8 == 0) {
            // Each bit in the bitmap corresponds to one byte in the bytemap
            // Note: the bitmap must be byte-aligned (8 bits)
            int bitmap = 0;
            for (unsigned int i = i_src; i < i_src + 8 && i < bytemap.size(); i++) {
                bitmap |= bytemap[i] << (i % 8);
            }
            bytemap[i_dst++] = bitmap;
        }
    }
}

void ColumnBuffer::to_bitmap(std::span<const uint8_t> bytemap, std::span<uint8_t> bitmap) {
    size_t i_dst = 0;
    for (size_t i_src = 0; i_src < bytemap.size(); ++i_src) {
        // Overwrite every 8 bytes with a one-byte bitmap
        if (i_src % 8 == 0) {
            // Each bit in the bitmap corresponds to one byte in the bytemap
            // Note: the bitmap must be byte-aligned (8 bits)
            bitmap[i_dst] = 0;
            for (size_t i = i_src; i < i_src + 8 && i < bytemap.size(); i++) {
                bitmap[i_dst] |= bytemap[i] << (i % 8);
            }
            ++i_dst;
        }
    }
}

MemoryMode ColumnBuffer::memory_mode(const Config& config) {
    if (config.contains(CONFIG_KEY_MEMORY_MODE)) {
        std::string value = config.get(CONFIG_KEY_MEMORY_MODE);
        if (value == "performance") {
            return MemoryMode::PERFORMANCE;
        } else if (value == "efficiency") {
            return MemoryMode::EFFICIENCY;
        } else {
            throw TileDBSOMAError(fmt::format("[ColumnBuffer] Unknown memory mode specified '{}'", value));
        }
    }

    return DEFAULT_MEMORY_MODE;
}

#pragma endregion

#pragma region public non-static

ColumnBuffer::ColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    size_t num_cells,
    size_t max_num_cells,
    size_t data_size,
    size_t max_data_size,
    bool is_var,
    bool is_nullable,
    std::optional<Enumeration> enumeration,
    bool is_ordered,
    MemoryMode mode)
    : num_cells_(num_cells)
    , data_size_(data_size)
    , max_num_cells_(max_num_cells)
    , max_data_size_(max_data_size)
    , mode_(mode)
    , name_(name)
    , type_(type)
    , type_size_(tiledb::impl::type_size(type))
    , is_var_(is_var)
    , is_nullable_(is_nullable)
    , enumeration_(enumeration)
    , is_ordered_(is_ordered) {
}

ColumnBuffer::ColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    size_t max_num_cells,
    size_t max_data_size,
    bool is_var,
    bool is_nullable,
    std::optional<Enumeration> enumeration,
    bool is_ordered,
    MemoryMode mode)
    : ColumnBuffer(name, type, 0, max_num_cells, 0, max_data_size, is_var, is_nullable, enumeration, is_ordered, mode) {
}

ColumnBuffer::~ColumnBuffer() {
}

void ColumnBuffer::attach(Query& query, std::optional<Subarray> subarray) {
    auto is_write = query.query_type() == TILEDB_WRITE;
    auto schema = query.array().schema();
    bool is_dense = schema.array_type() == TILEDB_DENSE;
    auto is_dim = schema.domain().has_dimension(name_);
    auto use_subarray = is_write && is_dense && is_dim;

    if (use_subarray && !subarray.has_value()) {
        throw TileDBSOMAError(
            "[ColumnBuffer::attach] Subarray must be provided to ColumnBuffer "
            "to attach to Query");
    }

    return use_subarray ? attach_subarray(*subarray) : attach_buffer(query);
}

std::vector<std::vector<std::byte>> ColumnBuffer::binaries() const {
    if (!is_var_) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer][binaries] Column '{}' is not var sized.", name_));
    }

    std::vector<std::vector<std::byte>> result;

    auto data_ptr = data().data();
    auto offsets_view = offsets();

    for (size_t i = 0; i < num_cells_; ++i) {
        size_t data_sz = offsets_view[i + 1] - offsets_view[i];

        result.emplace_back(data_sz);
        std::memcpy(result[0].data(), data_ptr + offsets_view[i], data_sz);
    }

    return result;
}

std::vector<std::string> ColumnBuffer::strings() const {
    if (!is_var_) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer][binaries] Column '{}' is not var sized.", name_));
    }

    std::vector<std::string> result;
    result.reserve(size());

    auto data_ptr = data<char>().data();
    auto offsets_view = offsets();

    for (size_t i = 0; i < size(); ++i) {
        result.emplace_back(data_ptr + offsets_view[i], offsets_view[i + 1] - offsets_view[i]);
    }

    return result;
}

std::string_view ColumnBuffer::string_view(size_t index) const {
    auto data_ptr = data<char>().data();
    auto offsets_view = offsets();

    return std::string_view(data_ptr + offsets_view[index], offsets_view[index + 1] - offsets_view[index]);
}

std::unique_ptr<IArrowBufferStorage> ColumnBuffer::export_buffers() {
    std::unique_ptr<uint8_t[]> validity_buffer = nullptr;

    if (is_nullable()) {
        size_t bitmap_size = (num_cells_ + 7) / 8;
        validity_buffer = std::make_unique_for_overwrite<uint8_t[]>(bitmap_size);

        ColumnBuffer::to_bitmap(this->validity(), std::span<uint8_t>(validity_buffer.get(), bitmap_size));
    }

    std::unique_ptr<std::byte[]> data_buffer = nullptr;
    if (type_ == TILEDB_BOOL) {
        size_t bitmap_size = (num_cells_ + 7) / 8;
        data_buffer = std::make_unique_for_overwrite<std::byte[]>(bitmap_size);

        ColumnBuffer::to_bitmap(
            std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(this->data().data()), num_cells_),
            std::span<uint8_t>(reinterpret_cast<uint8_t*>(data_buffer.get()), bitmap_size));
    } else {
        data_buffer = std::make_unique_for_overwrite<std::byte[]>(data_size_);
        std::memcpy(data_buffer.get(), this->data().data(), data_size_);
    }

    if (is_var()) {
        std::unique_ptr<uint64_t[]> offsets_buffer = std::make_unique_for_overwrite<uint64_t[]>(num_cells_ + 1);
        std::memcpy(offsets_buffer.get(), this->offsets().data(), (num_cells_ + 1) * sizeof(uint64_t));

        return std::make_unique<ArrayArrowBufferStorage>(
            std::move(data_buffer), std::move(offsets_buffer), num_cells_, std::move(validity_buffer));
    } else {
        return std::make_unique<ArrayArrowBufferStorage>(
            type(), num_cells_, std::move(data_buffer), std::move(validity_buffer));
    }
}

#pragma endregion

#pragma region private non-static

void ColumnBuffer::attach_buffer(Query& query) {
    query.set_data_buffer(name_, (void*)data().data(), max_data_size_ / type_size_);
    if (is_var_) {
        // Remove one offset for TileDB, which checks that the offsets and
        // validity buffers are the same size
        query.set_offsets_buffer(name_, const_cast<uint64_t*>(offsets().data()), max_num_cells_);
    }
    if (is_nullable_) {
        query.set_validity_buffer(name_, const_cast<uint8_t*>(validity().data()), max_num_cells_);
    }
}

void ColumnBuffer::attach_subarray(Subarray& subarray) const {
    switch (type_) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
        case TILEDB_BLOB:
            return attach_range<std::string>(subarray);
        case TILEDB_FLOAT32:
            return attach_range<float>(subarray);
        case TILEDB_FLOAT64:
            return attach_range<double>(subarray);
        case TILEDB_UINT8:
            return attach_range<uint8_t>(subarray);
        case TILEDB_INT8:
            return attach_range<int8_t>(subarray);
        case TILEDB_UINT16:
            return attach_range<uint16_t>(subarray);
        case TILEDB_INT16:
            return attach_range<int16_t>(subarray);
        case TILEDB_UINT32:
            return attach_range<uint32_t>(subarray);
        case TILEDB_INT32:
            return attach_range<int32_t>(subarray);
        case TILEDB_UINT64:
            return attach_range<uint64_t>(subarray);
        case TILEDB_INT64:
        case TILEDB_TIME_SEC:
        case TILEDB_TIME_MS:
        case TILEDB_TIME_US:
        case TILEDB_TIME_NS:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
            return attach_range<int64_t>(subarray);
        default:
            return;
    }
}

#pragma endregion

#pragma endregion

#pragma region ReadColumnBuffer

#pragma region public non-static

ReadColumnBuffer::~ReadColumnBuffer() {
}

size_t ReadColumnBuffer::update_size(const Query& query) {
    auto [num_offsets, num_elements] = query.result_buffer_elements()[std::string(name())];

    if (is_var()) {
        num_cells_ = num_offsets;
        // Set the extra offset value for arrow.
        offsets()[num_offsets] = data_size_ = num_elements;
    } else {
        num_cells_ = num_elements;
        data_size_ = num_elements * tiledb::impl::type_size(type());
    }

    return num_cells_;
}

#pragma endregion

#pragma endregion

#pragma region CArrayColumnBuffer

#pragma region public non-static

CArrayColumnBuffer::CArrayColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    size_t num_cells,
    size_t num_bytes,
    bool is_var,
    bool is_nullable,
    std::optional<Enumeration> enumeration,
    bool is_ordered,
    MemoryMode mode)
    : ReadColumnBuffer(name, type, num_cells, num_bytes, is_var, is_nullable, enumeration, is_ordered, mode) {
    LOG_DEBUG(fmt::format("[CArrayColumnBuffer] '{}' {} bytes", name, num_bytes));

    data_ = std::make_unique_for_overwrite<std::byte[]>(num_bytes);
    if (is_var) {
        offsets_ = std::make_unique_for_overwrite<uint64_t[]>(num_cells + 1);
    }
    if (is_nullable) {
        validity_ = std::make_unique_for_overwrite<uint8_t[]>(num_cells);
    }
}

CArrayColumnBuffer::~CArrayColumnBuffer() {
}

std::span<const std::byte> CArrayColumnBuffer::data() const {
    return std::span<const std::byte>(data_.get(), data_size_);
}

std::span<std::byte> CArrayColumnBuffer::data() {
    return std::span<std::byte>(data_.get(), data_size_);
}

std::span<const uint64_t> CArrayColumnBuffer::offsets() const {
    if (!is_var()) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer] Offsets buffer not defined for '{}'", name()));
    }

    return std::span<const uint64_t>(offsets_.get(), num_cells_ + 1);
}

std::span<uint64_t> CArrayColumnBuffer::offsets() {
    if (!is_var()) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer] Offsets buffer not defined for '{}'", name()));
    }

    return std::span<uint64_t>(offsets_.get(), num_cells_ + 1);
}

std::span<const uint8_t> CArrayColumnBuffer::validity() const {
    if (!is_nullable()) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer] Validity buffer not defined for '{}'", name()));
    }

    return std::span<const uint8_t>(validity_.get(), num_cells_);
}

std::span<uint8_t> CArrayColumnBuffer::validity() {
    if (!is_nullable()) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer] Validity buffer not defined for '{}'", name()));
    }

    return std::span<uint8_t>(validity_.get(), num_cells_);
}

std::unique_ptr<std::byte[]> CArrayColumnBuffer::release_and_reallocate_data() {
    std::unique_ptr<std::byte[]> data_buffer = std::make_unique_for_overwrite<std::byte[]>(max_data_size_);

    data_.swap(data_buffer);

    return data_buffer;
}

std::unique_ptr<uint64_t[]> CArrayColumnBuffer::release_and_reallocate_offsets() {
    std::unique_ptr<uint64_t[]> offset_buffer = std::make_unique_for_overwrite<uint64_t[]>(max_num_cells_);

    offsets_.swap(offset_buffer);

    return offset_buffer;
}

std::unique_ptr<IArrowBufferStorage> CArrayColumnBuffer::export_buffers() {
    // Validity bitmap is constructed the same way regardless of the `MemoryMode` selected
    std::unique_ptr<uint8_t[]> validity_buffer = nullptr;
    if (is_nullable()) {
        size_t bitmap_size = (num_cells_ + 7) / 8;
        validity_buffer = std::make_unique_for_overwrite<uint8_t[]>(bitmap_size);

        ColumnBuffer::to_bitmap(
            validity(), std::span<uint8_t>(reinterpret_cast<uint8_t*>(validity_buffer.get()), bitmap_size));
    }

    // Bool typed column need to be casted from bytemap to bitmap regardless of the `MemoryMode`
    std::unique_ptr<std::byte[]> data_buffer = nullptr;
    if (type() == TILEDB_BOOL) {
        size_t bitmap_size = (num_cells_ + 7) / 8;
        data_buffer = std::make_unique_for_overwrite<std::byte[]>(bitmap_size);

        ColumnBuffer::to_bitmap(
            std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(this->data().data()), num_cells_),
            std::span<uint8_t>(reinterpret_cast<uint8_t*>(data_buffer.get()), bitmap_size));
    }

    if (mode_ == MemoryMode::PERFORMANCE) {
        // Allocate a new buffer and swap with the current buffer containing the read data

        if (type() != TILEDB_BOOL) {
            data_buffer = std::make_unique_for_overwrite<std::byte[]>(max_data_size_);
            this->data_.swap(data_buffer);
        }

        if (is_var()) {
            std::unique_ptr<uint64_t[]> offsets_buffer = std::make_unique_for_overwrite<uint64_t[]>(max_num_cells_ + 1);
            this->offsets_.swap(offsets_buffer);

            return std::make_unique<ArrayArrowBufferStorage>(
                std::move(data_buffer), std::move(offsets_buffer), num_cells_, std::move(validity_buffer));
        } else {
            return std::make_unique<ArrayArrowBufferStorage>(
                type(), num_cells_, std::move(data_buffer), std::move(validity_buffer));
        }
    } else {
        if (type() != TILEDB_BOOL) {
            data_buffer = std::make_unique_for_overwrite<std::byte[]>(data_size_);

            if (data_size_ == max_data_size_) {
                // If the data buffer is filled completetly we move it to the arrow structure and reserve a new one for the `ColumnBuffer`
                this->data_.swap(data_buffer);
            } else {
                // Else copy the buffer contents to an appropriatly sized buffer to be used by Arrow
                std::memcpy(data_buffer.get(), this->data_.get(), data_size_);
            }
        }

        if (is_var()) {
            std::unique_ptr<uint64_t[]> offsets_buffer = std::make_unique_for_overwrite<uint64_t[]>(num_cells_ + 1);

            if (num_cells_ == max_num_cells_) {
                // If the offset buffer is filled completetly we move it to the arrow structure and reserve a new one for the `ColumnBuffer`
                this->offsets_.swap(offsets_buffer);
            } else {
                // Else copy the buffer contents to an appropriatly sized buffer to be used by Arrow
                std::memcpy(offsets_buffer.get(), this->offsets_.get(), (num_cells_ + 1) * sizeof(uint64_t));
            }

            return std::make_unique<ArrayArrowBufferStorage>(
                std::move(data_buffer), std::move(offsets_buffer), num_cells_, std::move(validity_buffer));
        } else {
            return std::make_unique<ArrayArrowBufferStorage>(
                type(), num_cells_, std::move(data_buffer), std::move(validity_buffer));
        }
    }
}

#pragma endregion

#pragma endregion

#pragma region ViewColumnBuffer

#pragma region public non-static

WriteColumnBuffer::~WriteColumnBuffer() {
}

std::span<const std::byte> WriteColumnBuffer::data() const {
    return std::span<const std::byte>(data_view_ != nullptr ? data_view_ : data_buffer_.get(), data_size_);
}

std::span<const uint64_t> WriteColumnBuffer::offsets() const {
    if (!is_var()) {
        throw TileDBSOMAError(fmt::format("[WriteColumnBuffer] Offsets buffer not defined for '{}'", name()));
    }

    return std::span<const uint64_t>(offsets_view_ != nullptr ? offsets_view_ : offsets_buffer_.get(), num_cells_ + 1);
}

std::span<const uint8_t> WriteColumnBuffer::validity() const {
    if (!is_nullable()) {
        throw TileDBSOMAError(fmt::format("[WriteColumnBuffer] Validity buffer not defined for '{}'", name()));
    }

    return std::span<const uint8_t>(validity_buffer_.get(), num_cells_);
}

#pragma endregion

#pragma endregion

//===================================================================
//= public static
//===================================================================

std::shared_ptr<ColumnBuffer> VectorColumnBuffer::create(std::shared_ptr<Array> array, std::string_view name) {
    auto schema = array->schema();
    auto name_str = std::string(name);  // string for TileDB API

    if (schema.has_attribute(name_str)) {
        auto attr = schema.attribute(name_str);
        auto type = attr.type();
        bool is_var = attr.cell_val_num() == TILEDB_VAR_NUM;
        bool is_nullable = attr.nullable();
        auto enum_name = AttributeExperimental::get_enumeration_name(schema.context(), attr);
        std::optional<Enumeration> enumeration = std::nullopt;
        bool is_ordered = false;
        if (enum_name.has_value()) {
            auto enmr = ArrayExperimental::get_enumeration(schema.context(), *array, *enum_name);
            is_ordered = enmr.ordered();
            enumeration = std::make_optional<Enumeration>(enmr);
        }

        if (!is_var && attr.cell_val_num() != 1) {
            throw TileDBSOMAError("[ColumnBuffer] Values per cell > 1 is not supported: " + name_str);
        }

        return VectorColumnBuffer::alloc(
            schema.context().config(), name_str, type, is_var, is_nullable, enumeration, is_ordered);

    } else if (schema.domain().has_dimension(name_str)) {
        auto dim = schema.domain().dimension(name_str);
        auto type = dim.type();
        bool is_var = dim.cell_val_num() == TILEDB_VAR_NUM || dim.type() == TILEDB_STRING_ASCII ||
                      dim.type() == TILEDB_STRING_UTF8;

        if (!is_var && dim.cell_val_num() != 1) {
            throw TileDBSOMAError("[ColumnBuffer] Values per cell > 1 is not supported: " + name_str);
        }

        return VectorColumnBuffer::alloc(schema.context().config(), name_str, type, is_var, false, std::nullopt, false);
    }

    throw TileDBSOMAError("[ColumnBuffer] Column name not found: " + name_str);
}

//===================================================================
//= public non-static
//===================================================================

VectorColumnBuffer::VectorColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    size_t num_cells,
    size_t num_bytes,
    bool is_var,
    bool is_nullable,
    std::optional<Enumeration> enumeration,
    bool is_ordered,
    MemoryMode mode)
    : ReadColumnBuffer(name, type, num_cells, num_bytes, is_var, is_nullable, enumeration, is_ordered, mode) {
    LOG_DEBUG(
        fmt::format(
            "[VectorColumnBuffer] '{}' {} bytes is_var={} is_nullable={}", name, num_bytes, is_var, is_nullable));
    // Call reserve() to allocate memory without initializing the contents.
    // This reduce the time to allocate the buffer and reduces the
    // resident memory footprint of the buffer.
    data_.resize(num_bytes);
    if (is_var) {
        offsets_.resize(num_cells + 1);  // extra offset for arrow
    }
    if (is_nullable) {
        validity_.resize(num_cells);
    }
}

VectorColumnBuffer::~VectorColumnBuffer() {
    LOG_TRACE(fmt::format("[ColumnBuffer] release '{}'", name()));
}

std::span<const std::byte> VectorColumnBuffer::data() const {
    return std::span<const std::byte>(data_.data(), data_size_);
}

std::span<std::byte> VectorColumnBuffer::data() {
    return std::span<std::byte>(data_.data(), data_size_);
}

std::span<const uint64_t> VectorColumnBuffer::offsets() const {
    if (!is_var()) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer] Offsets buffer not defined for '{}'", name()));
    }

    return std::span<const uint64_t>(offsets_.data(), num_cells_);
}

std::span<uint64_t> VectorColumnBuffer::offsets() {
    if (!is_var()) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer] Offsets buffer not defined for '{}'", name()));
    }

    return std::span<uint64_t>(offsets_.data(), num_cells_);
}

std::span<const uint8_t> VectorColumnBuffer::validity() const {
    if (!is_nullable()) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer] Validity buffer not defined for '{}'", name()));
    }

    return std::span<const uint8_t>(validity_.data(), num_cells_);
}

std::span<uint8_t> VectorColumnBuffer::validity() {
    if (!is_nullable()) {
        throw TileDBSOMAError(fmt::format("[ColumnBuffer] Validity buffer not defined for '{}'", name()));
    }

    return std::span<uint8_t>(validity_.data(), num_cells_);
}

std::unique_ptr<IArrowBufferStorage> VectorColumnBuffer::export_buffers() {
    // Validity bitmap is constructed the same way regardless of the `MemoryMode` selected
    std::vector<uint8_t, NoInitAlloc<uint8_t>> validity_buffer;
    if (is_nullable()) {
        size_t bitmap_size = (num_cells_ + 7) / 8;
        validity_buffer.resize(bitmap_size);

        ColumnBuffer::to_bitmap(
            validity(), std::span<uint8_t>(reinterpret_cast<uint8_t*>(validity_buffer.data()), bitmap_size));
    }

    // Bool typed column need to be casted from bytemap to bitmap regardless of the `MemoryMode`
    std::vector<std::byte, NoInitAlloc<std::byte>> data_buffer;
    if (type() == TILEDB_BOOL) {
        size_t bitmap_size = (num_cells_ + 7) / 8;
        data_buffer.resize(bitmap_size);

        ColumnBuffer::to_bitmap(
            std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(this->data().data()), num_cells_),
            std::span<uint8_t>(reinterpret_cast<uint8_t*>(data_buffer.data()), bitmap_size));
    }

    if (mode_ == MemoryMode::PERFORMANCE) {
        // Allocate a new buffer and swap with the current buffer containing the read data

        if (type() != TILEDB_BOOL) {
            data_buffer.resize(max_data_size_);
            this->data_.swap(data_buffer);
        }

        if (is_var()) {
            std::vector<uint64_t, NoInitAlloc<uint64_t>> offsets_buffer(max_num_cells_ + 1);
            this->offsets_.swap(offsets_buffer);

            return std::make_unique<VectorArrowBufferStorage>(
                std::move(data_buffer), std::move(offsets_buffer), num_cells_, std::move(validity_buffer));
        } else {
            return std::make_unique<VectorArrowBufferStorage>(
                type(), num_cells_, std::move(data_buffer), std::move(validity_buffer));
        }
    } else {
        if (type() != TILEDB_BOOL) {
            data_buffer.resize(data_size_);

            if (data_size_ == max_data_size_) {
                // If the data buffer is filled completetly we move it to the arrow structure and reserve a new one for the `ColumnBuffer`
                this->data_.swap(data_buffer);
            } else {
                // Else copy the buffer contents to an appropriatly sized buffer to be used by Arrow
                std::memcpy(data_buffer.data(), this->data_.data(), data_size_);
            }
        }

        if (is_var()) {
            std::vector<uint64_t, NoInitAlloc<uint64_t>> offsets_buffer(num_cells_ + 1);

            if (num_cells_ == max_num_cells_) {
                // If the offset buffer is filled completetly we move it to the arrow structure and reserve a new one for the `ColumnBuffer`
                this->offsets_.swap(offsets_buffer);
            } else {
                // Else copy the buffer contents to an appropriatly sized buffer to be used by Arrow
                std::memcpy(offsets_buffer.data(), this->offsets_.data(), (num_cells_ + 1) * sizeof(uint64_t));
            }

            return std::make_unique<VectorArrowBufferStorage>(
                std::move(data_buffer), std::move(offsets_buffer), num_cells_, std::move(validity_buffer));
        } else {
            return std::make_unique<VectorArrowBufferStorage>(
                type(), num_cells_, std::move(data_buffer), std::move(validity_buffer));
        }
    }
}

//===================================================================
//= private static
//===================================================================

std::shared_ptr<ColumnBuffer> VectorColumnBuffer::alloc(
    Config config,
    std::string_view name,
    tiledb_datatype_t type,
    bool is_var,
    bool is_nullable,
    std::optional<Enumeration> enumeration,
    bool is_ordered) {
    // Set number of bytes for the data buffer. Override with a value from
    // the config if present.
    auto num_bytes = DEFAULT_ALLOC_BYTES;
    if (config.contains(CONFIG_KEY_INIT_BYTES)) {
        auto value_str = config.get(CONFIG_KEY_INIT_BYTES);
        try {
            num_bytes = std::stoull(value_str);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(
                fmt::format("[ColumnBuffer] Error parsing {}: '{}' ({})", CONFIG_KEY_INIT_BYTES, value_str, e.what()));
        }
    }

    MemoryMode mode = ColumnBuffer::memory_mode(config);

    // bool is_dense = schema.array_type() == TILEDB_DENSE;
    // if (is_dense) {
    //     // TODO: Handle dense arrays similar to tiledb python module
    // }

    // For variable length column types, allocate an extra num_bytes to hold
    //   offset values. The number of cells is the set by the size of the
    //   offset type.
    // For non-variable length column types, the number of cells is computed
    //   from the type size.
    size_t num_cells = is_var ? num_bytes / sizeof(uint64_t) : num_bytes / tiledb::impl::type_size(type);

    return std::make_shared<VectorColumnBuffer>(
        name, type, num_cells, num_bytes, is_var, is_nullable, enumeration, is_ordered, mode);
}

void ColumnBuffer::resize(size_t num_bytes, bool preserve_data) {
    std::vector<std::byte> data_buffer(num_bytes);
    std::vector<uint64_t> offsets_buffer;
    std::vector<uint8_t> validity_buffer;

    size_t new_num_cells = is_var_ ? num_bytes / sizeof(uint64_t) : num_bytes / tiledb::impl::type_size(type());

    if (is_var_) {
        offsets_buffer.resize(new_num_cells + 1);
    }

    if (is_nullable_) {
        validity_buffer.resize(new_num_cells);
    }

    if (preserve_data) {
        size_t copy_num_bytes = is_var_ ? std::min(num_bytes, (num_cells_ != 0 ? offsets_[num_cells_] : 0)) :
                                          std::min(num_bytes, num_cells_ * impl::type_size(type()));

        std::memcpy(data_buffer.data(), data_.data(), copy_num_bytes);

        if (is_var_) {
            std::memcpy(
                offsets_buffer.data(), offsets_.data(), std::min(new_num_cells + 1, num_cells_ + 1) * sizeof(uint64_t));
        }

        if (is_nullable_) {
            std::memcpy(
                offsets_buffer.data(), offsets_.data(), std::min(new_num_cells + 1, num_cells_ + 1) * sizeof(uint64_t));
        }
    }

    data_ = data_buffer;
    offsets_ = offsets_buffer;
    validity_ = validity_buffer;

    num_cells_ = std::min(num_cells_, new_num_cells);
}

size_t ColumnBuffer::max_size() const {
    return data_.capacity();
}

size_t ColumnBuffer::max_num_cells() const {
    return is_var_ ? max_size() / sizeof(uint64_t) : max_size() / tiledb::impl::type_size(type());
}

}  // namespace tiledbsoma
