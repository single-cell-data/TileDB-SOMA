/**
 * @file   column_buffer.h
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

#ifndef COMMON_COLUMN_BUFFER_H
#define COMMON_COLUMN_BUFFER_H

#include <algorithm>
#include <concepts>
#include <span>
#include <stdexcept>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include <cmath>

#include "../arrow/interface.h"
#include "../common.h"
#include "../logging/logger.h"

namespace tiledbsoma::common {

enum class MemoryMode { PERFORMANCE, EFFICIENCY };

inline constexpr std::string_view CONFIG_KEY_MEMORY_MODE{"tiledb.builder.memory_mode"};

class ColumnBuffer {
   public:
    static MemoryMode memory_mode(const tiledb::Config& config);

    inline static const MemoryMode DEFAULT_MEMORY_MODE = MemoryMode::PERFORMANCE;

    ColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        uint64_t num_cells,
        uint64_t max_num_cells,
        uint64_t data_size,
        uint64_t max_data_size,
        bool is_var = false,
        bool is_nullable = false,
        std::optional<tiledb::Enumeration> enumeration = std::nullopt,
        MemoryMode mode = DEFAULT_MEMORY_MODE);

    ColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        uint64_t max_num_cells,
        uint64_t max_data_size,
        bool is_var = false,
        bool is_nullable = false,
        std::optional<tiledb::Enumeration> enumeration = std::nullopt,
        MemoryMode mode = DEFAULT_MEMORY_MODE);

    ColumnBuffer() = delete;
    ColumnBuffer(const ColumnBuffer&) = delete;
    ColumnBuffer(ColumnBuffer&&) = default;

    virtual ~ColumnBuffer();

    /**
     * @brief Attach the internal buffer the the given query. If a subarray is provided, attach the subarray insted of the internal buffers.
     * 
     * @param query The query to attach the buffers to
     * @param subarray Optional subarray to attach to query instead of the columns internal buffers. 
     *  This should only be provided for dense write queries and only for dimensions.
     */
    void attach(tiledb::Query& query, std::optional<tiledb::Subarray> subarray = std::nullopt);

    /**
     * @brief Return the number of cells taht are stored in the buffers.
     *
     * @return uint64_t
     */
    uint64_t cell_count() const;

    /**
     * @brief Return the maximum number of cells that can be stored in the buffers.
     *
     * @return uint64_t
     */
    uint64_t cell_capacity() const;

    /**
     * @brief Return the size of the data buffer in bytes.
     *
     * @return uint64_t
     */
    uint64_t data_size() const;

    /**
     * @brief Return the capacity of the data buffer in bytes.
     *
     * @return uint64_t
     */
    uint64_t data_capacity() const;

    /**
     * @brief Return a read-only view of the ColumnBuffer data.
     *
     * @tparam T Data type
     * @return std::span<const T> data view
     */
    template <typename T>
    std::span<const T> data() const {
        return std::span<const T>(reinterpret_cast<const T*>(data().data()), cell_count());
    }

    /**
     * @brief Return data as a vector of spans of bytes.
     *
     * @return std::vector<std::span<const std::byte>>
     */
    std::vector<std::span<const std::byte>> binaries() const;

    /**
     * @brief Return data in a vector of string_views.
     *
     * @return std::vector<std::string_view>
     */
    std::vector<std::string_view> strings() const;

    /**
     * @brief Return a string_view of the string at the provided cell index.
     *
     * @param index Cell index
     * @return std::string_view string view
     */
    std::string_view string_view(size_t index) const;

    /**
     * @brief Return a read-only view of the ColumnBuffer raw data.
     *
     * @return std::span<const std::byte> raw data view
     */
    virtual std::span<const std::byte> data() const = 0;

    /**
     * @brief Return a read-only view of the ColumnBuffer offsets.
     *
     * @return std::span<const uint64_t> offsets view
     */
    virtual std::span<const uint64_t> offsets() const = 0;

    /**
     * @brief Return a read-only view of the validity buffer.
     *
     * @return std::span<const uint8_t> validity view
     */
    virtual std::span<const uint8_t> validity() const = 0;

    /**
     * @brief Return the name of the buffer.
     *
     * @return std::string_view
     */
    std::string_view name() const;

    /**
     * @brief Return the type of the buffer.
     *
     * @return tiledb_datatype_t type
     */
    tiledb_datatype_t type() const;

    /**
     * @brief Return true if the buffer contains variable length data.
     */
    bool is_var() const;

    /**
     * @brief Return true if the buffer contains nullable data.
     */
    bool is_nullable() const;

    /**
     * @brief Get the enumeration associated with this column.
     */
    std::optional<tiledb::Enumeration> get_enumeration() const;

    /**
     * @brief Return true if the buffer contains an ordered enumeration.
     */
    bool is_ordered() const;

    /**
     * @brief Return the percentage of the data buffer that is filled.
     */
    double_t data_load_factor() const;

    /**
     * @brief The percentage of the offset and validity buffers that is filled if
     * any.
     */
    double_t cell_load_factor() const;

    /**
     * @brief Export the internal buffers as a struct able to be used to construct an Arrow array.
     */
    virtual std::unique_ptr<arrow::IArrowBufferStorage> export_buffers();

   protected:
    // Number of cells stored in this column
    uint64_t num_cells_;

    // Number of bytes used from the data buffer
    uint64_t data_size_;

    // Maximum number of cells this column can hold
    uint64_t max_num_cells_;

    // Maximum number of bytes the data buffer can hold
    uint64_t max_data_size_;

    // Dictates how to handle internal buffers when converting from `ColumnBuffer` to `IArrowBufferStorage`
    MemoryMode mode_;

   private:
    /**
     * @brief Attach the internal buffers to the query object
     * 
     * @param query The query to attach to.
     */
    void attach_buffer(tiledb::Query& query);

    /**
     * @brief Set the column value range to the subarray object.
     * 
     * @param subarray The subarray to set the range to.
     */
    void attach_subarray(tiledb::Subarray& subarray) const;

    // Name of the column from the schema.
    std::string name_;

    // Data type of the column from the schema.
    tiledb_datatype_t type_;

    // Bytes per element.
    uint64_t type_size_;

    // If true, the data type is variable length
    bool is_var_;

    // If true, the data is nullable
    bool is_nullable_;

    // If applicable, the Enumeration associated with the column
    std::optional<tiledb::Enumeration> enumeration_;
};

class ReadColumnBuffer : public ColumnBuffer {
   public:
    using ColumnBuffer::ColumnBuffer;

    ReadColumnBuffer() = delete;
    ReadColumnBuffer(const ReadColumnBuffer&) = delete;
    ReadColumnBuffer(ReadColumnBuffer&&) = default;

    virtual ~ReadColumnBuffer();

    /**
   * @brief Size num_cells_ to match the read query results.
   *
   * @param query TileDB query
   */
    uint64_t update_size(const tiledb::Query& query);

    /**
     * @brief
     * 
     * @param num_bytes
     * @param num_cells
     * @param retain_data
     */
    virtual void resize(uint64_t num_bytes, uint64_t num_cells, bool retain_data = false) = 0;

    /**
   * @brief Return a view of the ColumnBuffer data.
   *
   * @tparam T Data type
   * @return std::span<T> data view
   */
    template <typename T>
    std::span<T> data() {
        return std::span<T>(reinterpret_cast<T*>(data().data()), cell_count());
    }

    /**
   * @brief Return a view of the ColumnBuffer raw data.
   *
   * @return std::span<const std::byte> raw data view
   */
    virtual std::span<std::byte> data() = 0;

    /**
   * @brief Return a view of the ColumnBuffer offsets.
   *
   * @return std::span<const uint64_t> offsets view
   */
    virtual std::span<uint64_t> offsets() = 0;

    /**
   * @brief Return a view of the validity buffer.
   *
   * @return std::span<const uint8_t> validity view
   */
    virtual std::span<uint8_t> validity() = 0;

    using ColumnBuffer::data;
    using ColumnBuffer::offsets;
    using ColumnBuffer::validity;
};

class CArrayColumnBuffer : public ReadColumnBuffer {
   public:
    /**
   * @brief Construct a new ColumnBuffer object
   *
   * @param name Column name
   * @param type TileDB datatype
   * @param num_cells Number of cells to allocate for offsets and validity
   * @param num_bytes Number of bytes to allocate for data
   * @param is_var Column type is variable length
   * @param is_nullable Column can contain null values
   * @param enumeration Optional Enumeration associated with column
   * @param is_ordered Optional Enumeration is ordered
   */
    CArrayColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        uint64_t num_cells,
        uint64_t num_bytes,
        bool is_var = false,
        bool is_nullable = false,
        std::optional<tiledb::Enumeration> enumeration = std::nullopt,
        MemoryMode mode = DEFAULT_MEMORY_MODE);

    CArrayColumnBuffer() = delete;
    CArrayColumnBuffer(const CArrayColumnBuffer&) = delete;
    CArrayColumnBuffer(CArrayColumnBuffer&&) = default;

    ~CArrayColumnBuffer();

    using ReadColumnBuffer::data;

    std::span<std::byte> data() override;
    std::span<uint64_t> offsets() override;
    std::span<uint8_t> validity() override;

    std::span<const std::byte> data() const override;
    std::span<const uint64_t> offsets() const override;
    std::span<const uint8_t> validity() const override;

    void resize(uint64_t num_bytes, uint64_t num_cells, bool retain_data = false) override;

    std::unique_ptr<std::byte[]> release_data();
    std::unique_ptr<uint64_t[]> release_offsets();

    std::unique_ptr<arrow::IArrowBufferStorage> export_buffers() override;

   private:
    std::unique_ptr<std::byte[]> data_;
    std::unique_ptr<uint64_t[]> offsets_;
    std::unique_ptr<uint8_t[]> validity_;
};

class WriteColumnBuffer : public ColumnBuffer {
   public:
    /**
   * @brief Create a ColumnBuffer from an array and column name.
   *
   * @param array TileDB array
   * @param name TileDB dimension or attribute name
   * @return IColumnBuffer
   */
    template <typename DataStorage, typename OffsetStorage>
        requires is_data_buffer<DataStorage> && is_offset_buffer<OffsetStorage>
    static std::shared_ptr<ColumnBuffer> create(
        const tiledb::Array& array,
        std::string_view name,
        uint64_t num_elements,
        DataStorage data_buffer,
        OffsetStorage offsets_buffer,
        std::unique_ptr<uint8_t[]> validity_buffer,
        bool copy_buffers = false) {
        if (array.query_type() != TILEDB_WRITE) {
            throw std::runtime_error("[WriteColumnBuffer] WriteColumnBuffer should only be used during write.");
        }

        auto schema = array.schema();
        auto name_str = std::string(name);  // string for TileDB API

        if (schema.has_attribute(name_str)) {
            auto attr = schema.attribute(name_str);
            auto type = attr.type();
            bool is_var = attr.cell_val_num() == TILEDB_VAR_NUM;
            bool is_nullable = attr.nullable();

            if (!is_var && attr.cell_val_num() != 1) {
                throw std::runtime_error("[WriteColumnBuffer] Values per cell > 1 is not supported: " + name_str);
            }

            uint64_t num_bytes = is_var ? offsets_buffer[num_elements] : num_elements * tiledb::impl::type_size(type);

            return std::make_shared<WriteColumnBuffer>(
                name,
                type,
                num_elements,
                num_bytes,
                is_var,
                is_nullable,
                std::move(data_buffer),
                std::move(offsets_buffer),
                std::move(validity_buffer),
                copy_buffers);

        } else if (schema.domain().has_dimension(name_str)) {
            auto dim = schema.domain().dimension(name_str);
            auto type = dim.type();
            bool is_var = dim.cell_val_num() == TILEDB_VAR_NUM || dim.type() == TILEDB_STRING_ASCII ||
                          dim.type() == TILEDB_STRING_UTF8;

            if (!is_var && dim.cell_val_num() != 1) {
                throw std::runtime_error("[WriteColumnBuffer] Values per cell > 1 is not supported: " + name_str);
            }

            uint64_t num_bytes = is_var ? offsets_buffer[num_elements] : num_elements * tiledb::impl::type_size(type);

            return std::make_shared<WriteColumnBuffer>(
                name,
                type,
                num_elements,
                num_bytes,
                is_var,
                false,
                std::move(data_buffer),
                std::move(offsets_buffer),
                nullptr,
                copy_buffers);
        }

        throw std::runtime_error("[WriteColumnBuffer] Column name not found: " + name_str);
    }

    template <typename DataStorage, typename OffsetStorage>
        requires is_data_buffer<DataStorage> && is_offset_buffer<OffsetStorage>
    WriteColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        uint64_t num_cells,
        uint64_t num_bytes,
        bool is_var,
        bool is_nullable,
        DataStorage data_buffer,
        OffsetStorage offsets_buffer,
        std::unique_ptr<uint8_t[]> validity_buffer,
        bool copy_buffers = false)
        : ColumnBuffer(name, type, num_cells, num_cells, num_bytes, num_bytes, is_var, is_nullable, std::nullopt) {
        if (data_buffer == nullptr) {
            throw std::runtime_error(
                "[WriteColumnBuffer] Supplied data buffer is null for column '" + std::string(name) + "'");
        }

        if constexpr (std::is_same_v<std::unique_ptr<std::byte[]>, DataStorage>) {
            // Data buffer ownership is passed to the write column buffer
            data_buffer_ = std::move(data_buffer);
        } else {
            if (copy_buffers) {
                data_buffer_ = std::make_unique_for_overwrite<std::byte[]>(data_size_);
                std::memcpy(data_buffer_.get(), data_buffer, data_size_);
            } else {
                data_view_ = reinterpret_cast<const std::byte*>(data_buffer);
            }
        }

        if (is_var) {
            if (offsets_buffer == nullptr) {
                throw std::runtime_error(
                    "[WriteColumnBuffer] Supplied offset buffer is null for var sized column '" + std::string(name) +
                    "'");
            }

            if constexpr (std::is_same_v<std::unique_ptr<uint64_t[]>, OffsetStorage>) {
                // Offset buffer ownership is passed to the write column buffer
                offsets_buffer_ = std::move(offsets_buffer);
            } else if constexpr (std::is_same_v<uint32_t, std::remove_const_t<std::remove_pointer_t<OffsetStorage>>>) {
                // A new offset buffer will be contructed using the input offset buffer
                offsets_buffer_ = std::make_unique_for_overwrite<uint64_t[]>(num_cells_ + 1);

                // Fill the buffer with offset and cast to the correct type
                std::transform(
                    offsets_buffer, offsets_buffer + num_cells_ + 1, offsets_buffer_.get(), [](const uint32_t& offset) {
                        return static_cast<uint64_t>(offset);
                    });
            } else {
                if (copy_buffers) {
                    offsets_buffer_ = std::make_unique_for_overwrite<uint64_t[]>(num_cells_ + 1);
                    std::memcpy(offsets_buffer_.get(), offsets_buffer, (num_cells_ + 1) * sizeof(uint64_t));
                } else {
                    offsets_view_ = offsets_buffer;
                }
            }
        }

        if (is_nullable) {
            if (validity_buffer) {
                validity_buffer_ = std::move(validity_buffer);
            } else {
                validity_buffer_ = std::make_unique_for_overwrite<uint8_t[]>(num_cells_);
                std::fill(validity_buffer_.get(), validity_buffer_.get() + num_cells_, 1);
            }
        } else {
            // As of version 1.15.6 we were throwing here. However, we found a
            // compatibility issue with pyarrow versions below 17. Thus we log and
            // continue.
            if (validity_buffer) {
                logging::LOG_DEBUG(
                    "[WriteColumnBuffer] Validity buffer passed for column '" + std::string(name) +
                    "' is being ignored");
            }
        }
    }

    WriteColumnBuffer() = delete;
    WriteColumnBuffer(const WriteColumnBuffer&) = delete;
    WriteColumnBuffer(WriteColumnBuffer&&) = default;

    virtual ~WriteColumnBuffer();

    std::span<const std::byte> data() const override;
    std::span<const uint64_t> offsets() const override;
    std::span<const uint8_t> validity() const override;

   private:
    const std::byte* data_view_ = nullptr;
    const uint64_t* offsets_view_ = nullptr;

    std::unique_ptr<std::byte[]> data_buffer_;
    std::unique_ptr<uint64_t[]> offsets_buffer_;
    std::unique_ptr<uint8_t[]> validity_buffer_;
};

}  // namespace tiledbsoma::common
#endif
