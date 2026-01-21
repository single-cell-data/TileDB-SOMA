/**
 * @file   column_buffer.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc.
 *
 * @section DESCRIPTION
 *
 * This file defines the a ColumBuffer class.
 */

#include "column_buffer.h"
#include "../arrow/arrow_buffer.h"
#include "../arrow/utils.h"
#include "../logging/impl/logger.h"

namespace tiledbsoma::common {

#pragma region ColumnBuffer

#pragma region public static

MemoryMode ColumnBuffer::memory_mode(const tiledb::Config& config) {
    if (config.contains(CONFIG_KEY_MEMORY_MODE)) {
        std::string value = config.get(CONFIG_KEY_MEMORY_MODE.data());
        if (value == "performance") {
            return MemoryMode::PERFORMANCE;
        } else if (value == "efficiency") {
            return MemoryMode::EFFICIENCY;
        } else {
            logging::LOG_WARN(
                fmt::format(
                    "[ColumnBuffer] Unknown memory mode specified '{}'. Default mode "
                    "selected.",
                    value));
        }
    }

    return DEFAULT_MEMORY_MODE;
}

#pragma endregion

#pragma region public non-static

ColumnBuffer::ColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    uint64_t num_cells,
    uint64_t max_num_cells,
    uint64_t data_size,
    uint64_t max_data_size,
    bool is_var,
    bool is_nullable,
    std::optional<tiledb::Enumeration> enumeration,
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
    , enumeration_(enumeration) {
}

ColumnBuffer::ColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    uint64_t max_num_cells,
    uint64_t max_data_size,
    bool is_var,
    bool is_nullable,
    std::optional<tiledb::Enumeration> enumeration,
    MemoryMode mode)
    : ColumnBuffer(name, type, 0, max_num_cells, 0, max_data_size, is_var, is_nullable, enumeration, mode) {
}

ColumnBuffer::~ColumnBuffer() {
}

uint64_t ColumnBuffer::cell_count() const {
    return num_cells_;
}

uint64_t ColumnBuffer::cell_capacity() const {
    return max_num_cells_;
}

uint64_t ColumnBuffer::data_size() const {
    return data_size_;
}

uint64_t ColumnBuffer::data_capacity() const {
    return max_data_size_;
}

std::string_view ColumnBuffer::name() const {
    return name_;
}

tiledb_datatype_t ColumnBuffer::type() const {
    return type_;
}

bool ColumnBuffer::is_var() const {
    return is_var_;
}

bool ColumnBuffer::is_nullable() const {
    return is_nullable_;
}

std::optional<tiledb::Enumeration> ColumnBuffer::get_enumeration() const {
    return enumeration_;
}

bool ColumnBuffer::is_ordered() const {
    return enumeration_->ordered();
}

double_t ColumnBuffer::data_load_factor() const {
    return static_cast<double_t>(data_size_) / max_data_size_;
}

double_t ColumnBuffer::cell_load_factor() const {
    return static_cast<double_t>(num_cells_) / max_num_cells_;
}

void ColumnBuffer::attach(tiledb::Query& query, std::optional<tiledb::Subarray> subarray) {
    bool is_write = query.query_type() == TILEDB_WRITE;
    const tiledb::ArraySchema schema = query.array().schema();
    bool is_dense = schema.array_type() == TILEDB_DENSE;
    bool is_dim = schema.domain().has_dimension(name_);
    bool use_subarray = is_write && is_dense && is_dim;

    if (use_subarray && !subarray.has_value()) {
        throw std::runtime_error(
            "[ColumnBuffer::attach] Subarray must be provided to ColumnBuffer "
            "to attach to Query");
    }

    return use_subarray ? attach_subarray(*subarray) : attach_buffer(query);
}

std::vector<std::span<const std::byte>> ColumnBuffer::binaries() const {
    std::vector<std::span<const std::byte>> result;

    auto data_ptr = data().data();
    auto offsets_view = offsets();

    for (size_t i = 0; i < num_cells_; ++i) {
        result.emplace_back(data_ptr + offsets_view[i], offsets_view[i + 1] - offsets_view[i]);
    }

    return result;
}

std::vector<std::string_view> ColumnBuffer::strings() const {
    std::vector<std::string_view> result;
    result.reserve(cell_count());

    auto data_ptr = data<char>().data();
    auto offsets_view = offsets();

    for (size_t i = 0; i < cell_count(); ++i) {
        result.emplace_back(data_ptr + offsets_view[i], offsets_view[i + 1] - offsets_view[i]);
    }

    return result;
}

std::string_view ColumnBuffer::string_view(size_t index) const {
    auto data_ptr = data<char>().data();
    auto offsets_view = offsets();

    return std::string_view(data_ptr + offsets_view[index], offsets_view[index + 1] - offsets_view[index]);
}

std::unique_ptr<arrow::IArrowBufferStorage> ColumnBuffer::export_buffers() {
    std::unique_ptr<uint8_t[]> validity_buffer = nullptr;

    if (is_nullable()) {
        validity_buffer = arrow::bytemap_to_bitmap(validity().data(), cell_count(), 0);
    }

    std::unique_ptr<std::byte[]> data_buffer = nullptr;
    if (type_ == TILEDB_BOOL) {
        data_buffer = std::unique_ptr<std::byte[]>(
            reinterpret_cast<std::byte*>(arrow::bytemap_to_bitmap(data<uint8_t>().data(), cell_count(), 0).release()));
    } else {
        data_buffer = std::make_unique_for_overwrite<std::byte[]>(data_size());
        std::memcpy(data_buffer.get(), this->data().data(), data_size());
    }

    if (is_var()) {
        std::unique_ptr<uint64_t[]> offsets_buffer = std::make_unique_for_overwrite<uint64_t[]>(num_cells_ + 1);
        std::memcpy(offsets_buffer.get(), this->offsets().data(), (cell_count() + 1) * sizeof(uint64_t));

        return std::make_unique<arrow::ArrayArrowBufferStorage>(
            std::move(data_buffer), std::move(offsets_buffer), cell_count(), std::move(validity_buffer));
    } else {
        return std::make_unique<arrow::ArrayArrowBufferStorage>(
            type(), cell_count(), std::move(data_buffer), std::move(validity_buffer));
    }
}

#pragma endregion

#pragma region private non-static

void ColumnBuffer::attach_buffer(tiledb::Query& query) {
    query.set_data_buffer(name_, (void*)data().data(), data_capacity() / type_size_);
    if (is_var_) {
        query.set_offsets_buffer(name_, const_cast<uint64_t*>(offsets().data()), cell_capacity());
    }
    if (is_nullable_) {
        query.set_validity_buffer(name_, const_cast<uint8_t*>(validity().data()), cell_capacity());
    }
}

void ColumnBuffer::attach_subarray(tiledb::Subarray& subarray) const {
    auto attach_range = [&]<typename T>() {
        if constexpr (std::is_same_v<T, std::string>) {
            subarray.add_range(name_, std::string(string_view(0)), std::string(string_view(1)));
        } else {
            subarray.add_range(name_, data<T>()[0], data<T>()[1]);
        }
    };

    switch (type_) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
        case TILEDB_BLOB:
            return attach_range.template operator()<std::string>();
        case TILEDB_FLOAT32:
            return attach_range.template operator()<float>();
        case TILEDB_FLOAT64:
            return attach_range.template operator()<double>();
        case TILEDB_UINT8:
            return attach_range.template operator()<uint8_t>();
        case TILEDB_INT8:
            return attach_range.template operator()<int8_t>();
        case TILEDB_UINT16:
            return attach_range.template operator()<uint16_t>();
        case TILEDB_INT16:
            return attach_range.template operator()<int16_t>();
        case TILEDB_UINT32:
            return attach_range.template operator()<uint32_t>();
        case TILEDB_INT32:
            return attach_range.template operator()<int32_t>();
        case TILEDB_UINT64:
            return attach_range.template operator()<uint64_t>();
        case TILEDB_INT64:
        case TILEDB_TIME_SEC:
        case TILEDB_TIME_MS:
        case TILEDB_TIME_US:
        case TILEDB_TIME_NS:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
            return attach_range.template operator()<int64_t>();
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

uint64_t ReadColumnBuffer::update_size(const tiledb::Query& query) {
    auto [num_offsets, num_elements] = query.result_buffer_elements()[name().data()];

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
    uint64_t num_cells,
    uint64_t num_bytes,
    bool is_var,
    bool is_nullable,
    std::optional<tiledb::Enumeration> enumeration,
    MemoryMode mode)
    : ReadColumnBuffer(name, type, num_cells, num_bytes, is_var, is_nullable, enumeration, mode) {
    logging::LOG_DEBUG(fmt::format("[CArrayColumnBuffer] '{}' {} bytes", name, num_bytes));

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
    return std::span<const std::byte>(data_.get(), data_size());
}

std::span<std::byte> CArrayColumnBuffer::data() {
    return std::span<std::byte>(data_.get(), data_size());
}

std::span<const uint64_t> CArrayColumnBuffer::offsets() const {
    if (!is_var()) {
        throw std::runtime_error(fmt::format("[ColumnBuffer] Offsets buffer not defined for '{}'", name()));
    }

    return std::span<const uint64_t>(offsets_.get(), cell_count() + 1);
}

std::span<uint64_t> CArrayColumnBuffer::offsets() {
    if (!is_var()) {
        throw std::runtime_error(fmt::format("[ColumnBuffer] Offsets buffer not defined for '{}'", name()));
    }

    return std::span<uint64_t>(offsets_.get(), cell_count() + 1);
}

std::span<const uint8_t> CArrayColumnBuffer::validity() const {
    if (!is_nullable()) {
        throw std::runtime_error(fmt::format("[ColumnBuffer] Validity buffer not defined for '{}'", name()));
    }

    return std::span<const uint8_t>(validity_.get(), cell_count());
}

std::span<uint8_t> CArrayColumnBuffer::validity() {
    if (!is_nullable()) {
        throw std::runtime_error(fmt::format("[ColumnBuffer] Validity buffer not defined for '{}'", name()));
    }

    return std::span<uint8_t>(validity_.get(), cell_count());
}

void CArrayColumnBuffer::resize(uint64_t num_bytes, uint64_t num_cells, bool retain_data) {
    std::unique_ptr<std::byte[]> data_buffer = std::make_unique_for_overwrite<std::byte[]>(num_bytes);
    std::unique_ptr<uint64_t[]> offsets_buffer;
    std::unique_ptr<uint8_t[]> validity_buffer;

    if (is_var()) {
        offsets_buffer = std::make_unique_for_overwrite<uint64_t[]>(num_cells + 1);
    }

    if (is_nullable()) {
        validity_buffer = std::make_unique_for_overwrite<uint8_t[]>(num_cells);
    }

    if (retain_data) {
        std::memcpy(data_buffer.get(), data_.get(), std::min(num_bytes, data_size()));

        if (is_var()) {
            std::memcpy(offsets_buffer.get(), offsets_.get(), std::min(num_cells + 1, cell_count() + 1));
        }

        if (is_nullable()) {
            std::memcpy(validity_buffer.get(), validity_.get(), std::min(num_cells, cell_count()));
        }

        data_size_ = std::min(num_bytes, data_size());
        num_cells_ = std::min(num_cells, cell_count());
    } else {
        data_size_ = 0;
        num_cells_ = 0;
    }

    max_data_size_ = num_bytes;
    max_num_cells_ = num_cells;

    data_ = std::move(data_buffer);
    offsets_ = std::move(offsets_buffer);
    validity_ = std::move(validity_buffer);
}

std::unique_ptr<std::byte[]> CArrayColumnBuffer::release_data() {
    std::unique_ptr<std::byte[]> data_buffer = std::make_unique_for_overwrite<std::byte[]>(data_capacity());

    data_.swap(data_buffer);

    return data_buffer;
}

std::unique_ptr<uint64_t[]> CArrayColumnBuffer::release_offsets() {
    std::unique_ptr<uint64_t[]> offset_buffer = std::make_unique_for_overwrite<uint64_t[]>(cell_capacity());

    offsets_.swap(offset_buffer);

    return offset_buffer;
}

std::unique_ptr<arrow::IArrowBufferStorage> CArrayColumnBuffer::export_buffers() {
    // Validity bitmap is constructed the same way regardless of the `MemoryMode`
    // selected
    std::unique_ptr<uint8_t[]> validity_buffer = nullptr;
    if (is_nullable()) {
        validity_buffer = arrow::bytemap_to_bitmap(validity().data(), cell_count(), 0);
    }

    // Bool typed column need to be casted from bytemap to bitmap regardless of
    // the `MemoryMode`
    std::unique_ptr<std::byte[]> data_buffer = nullptr;
    if (type() == TILEDB_BOOL) {
        data_buffer = std::unique_ptr<std::byte[]>(
            reinterpret_cast<std::byte*>(arrow::bytemap_to_bitmap(data<uint8_t>().data(), cell_count(), 0).release()));
    }

    if (mode_ == MemoryMode::PERFORMANCE) {
        // Allocate a new buffer and swap with the current buffer containing the
        // read data

        if (type() != TILEDB_BOOL) {
            data_buffer = std::make_unique_for_overwrite<std::byte[]>(data_capacity());
            this->data_.swap(data_buffer);
        }

        if (is_var()) {
            std::unique_ptr<uint64_t[]> offsets_buffer = std::make_unique_for_overwrite<uint64_t[]>(
                cell_capacity() + 1);
            this->offsets_.swap(offsets_buffer);

            return std::make_unique<arrow::ArrayArrowBufferStorage>(
                std::move(data_buffer), std::move(offsets_buffer), cell_count(), std::move(validity_buffer));
        } else {
            return std::make_unique<arrow::ArrayArrowBufferStorage>(
                type(), cell_count(), std::move(data_buffer), std::move(validity_buffer));
        }
    } else {
        if (type() != TILEDB_BOOL) {
            data_buffer = std::make_unique_for_overwrite<std::byte[]>(data_size());

            if (data_size() == data_capacity()) {
                // If the data buffer is filled completetly we move it to the arrow
                // structure and reserve a new one for the `ColumnBuffer`
                this->data_.swap(data_buffer);
            } else {
                // Else copy the buffer contents to an appropriatly sized buffer to be
                // used by Arrow
                std::memcpy(data_buffer.get(), this->data_.get(), data_size());
            }
        }

        if (is_var()) {
            std::unique_ptr<uint64_t[]> offsets_buffer = std::make_unique_for_overwrite<uint64_t[]>(cell_count() + 1);

            if (cell_count() == cell_capacity()) {
                // If the offset buffer is filled completetly we move it to the arrow
                // structure and reserve a new one for the `ColumnBuffer`
                this->offsets_.swap(offsets_buffer);
            } else {
                // Else copy the buffer contents to an appropriatly sized buffer to be
                // used by Arrow
                std::memcpy(offsets_buffer.get(), this->offsets_.get(), (cell_count() + 1) * sizeof(uint64_t));
            }

            return std::make_unique<arrow::ArrayArrowBufferStorage>(
                std::move(data_buffer), std::move(offsets_buffer), cell_count(), std::move(validity_buffer));
        } else {
            return std::make_unique<arrow::ArrayArrowBufferStorage>(
                type(), cell_count(), std::move(data_buffer), std::move(validity_buffer));
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
    return data_;
}

std::span<const uint64_t> WriteColumnBuffer::offsets() const {
    if (!is_var()) {
        throw std::runtime_error(fmt::format("[WriteColumnBuffer] Offsets buffer not defined for '{}'", name()));
    }

    return offsets_;
}

std::span<const uint8_t> WriteColumnBuffer::validity() const {
    if (!is_nullable()) {
        throw std::runtime_error(fmt::format("[WriteColumnBuffer] Validity buffer not defined for '{}'", name()));
    }

    return validity_;
}

#pragma endregion

#pragma endregion

}  // namespace tiledbsoma::common
