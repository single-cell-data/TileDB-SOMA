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
 *   This declares the column buffer API
 */

#ifndef COLUMN_BUFFER_H
#define COLUMN_BUFFER_H

#include <algorithm>
#include <concepts>
#include <span>
#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include "../utils/arrow_adapter.h"
#include "../utils/common.h"
#include "logger_public.h"
#include "soma_context.h"

namespace tiledbsoma {

using namespace tiledb;

class ColumnBuffer {
   public:
    /**
     * @brief Convert a bytemap to a bitmap in place.
     *
     */
    static void to_bitmap(std::span<uint8_t> bytemap);

    /**
     * @brief Convert a bytemap to a bitmap.
     *
     */
    static void to_bitmap(std::span<const uint8_t> bytemap, std::span<uint8_t> bitmap);

    ColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        size_t max_num_cells,
        size_t data_size,
        size_t max_data_size,
        bool is_var = false,
        bool is_nullable = false,
        std::optional<Enumeration> enumeration = std::nullopt,
        bool is_ordered = false);

    ColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        size_t max_num_cells,
        size_t max_data_size,
        bool is_var = false,
        bool is_nullable = false,
        std::optional<Enumeration> enumeration = std::nullopt,
        bool is_ordered = false);

    ColumnBuffer() = delete;
    ColumnBuffer(const ColumnBuffer&) = delete;
    ColumnBuffer(ColumnBuffer&&) = default;

    ~ColumnBuffer();

    void attach(Query& query, std::optional<Subarray> subarray = std::nullopt);

    /**
     * @brief Return the number of cells in the buffer.
     *
     * @return size_t
     */
    size_t size() const {
        return num_cells_;
    }

    /**
     * @brief Return size of the data buffer.
     *
     * @return uint64_t
     */
    uint64_t data_size() const {
        return data_size_;
    }

    /**
     * @brief Return a read-only view of the ColumnBuffer data.
     *
     * @tparam T Data type
     * @return std::span<const T> data view
     */
    template <typename T>
    std::span<const T> data() const {
        return std::span<const T>(reinterpret_cast<const T*>(data().data()), size());
    }

    /**
     * @brief Return data in a vector of binary buffers.
     *
     * @return std::vector<std::vector<std::byte>>
     */
    std::vector<std::vector<std::byte>> binaries() const;

    /**
     * @brief Return data in a vector of strings.
     *
     * @return std::vector<std::string>
     */
    std::vector<std::string> strings() const;

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
    std::string_view name() const {
        return name_;
    }

    /**
     * @brief Return the type of the buffer.
     *
     * @return tiledb_datatype_t type
     */
    tiledb_datatype_t type() const {
        return type_;
    }

    /**
     * @brief Return true if the buffer contains variable length data.
     */
    bool is_var() const {
        return is_var_;
    }

    /**
     * @brief Return true if the buffer contains nullable data.
     */
    bool is_nullable() const {
        return is_nullable_;
    }

    std::optional<Enumeration> get_enumeration_info() const {
        return enumeration_;
    }

    /**
     * @brief Add an optional enumeration vector,
     *
     */
    void add_enumeration(const std::vector<std::string>& vec) {
        enums_ = vec;
        has_enumeration_ = true;
    }

    /**
     * @brief Return true if the buffer contains enumeration.
     */
    bool has_enumeration() const {
        return has_enumeration_;
    }

    /**
     * @brief Return optional enumeration vector,
     *
     */
    std::vector<std::string> get_enumeration() {
        return enums_;
    }

    /**
     * @brief Convert enumeration (aka dictionary)
     *
     */
    void convert_enumeration() {
        if (!has_enumeration_) {
            throw TileDBSOMAError("[ColumnBuffer] No enumeration defined for " + name_);
        }
        const size_t n_vec = enums_.size();  // plus one for extra offset
        enum_offsets_.resize(n_vec + 1);
        enum_str_ = "";
        uint32_t cumlen = 0;
        for (size_t i = 0; i < n_vec; i++) {
            std::string s(enums_[i]);
            enum_str_ += s;
            enum_offsets_[i] = cumlen;
            cumlen += static_cast<uint32_t>(s.length());
        }
        enum_offsets_[n_vec] = cumlen;
    }

    /**
     * @brief Return optional enumeration offsets vector
     *
     */
    std::span<uint32_t> enum_offsets() {
        if (!has_enumeration_) {
            throw TileDBSOMAError("[ColumnBuffer] No enumeration defined for " + name_);
        }
        return std::span<uint32_t>(enum_offsets_.data(), enum_offsets_.size());
    }

    /**
     * @brief Return optional enumeration string
     *
     */
    std::span<char> enum_string() {
        if (!has_enumeration_) {
            throw TileDBSOMAError("[ColumnBuffer] No enumeration defined for " + name_);
        }
        return std::span<char>(enum_str_.data(), enum_str_.length());
    }

    /**
     * @brief Return true if the buffer contains an ordered enumeration.
     */
    bool is_ordered() const {
        return is_ordered_;
    }

    /**
     * @brief Return the percentage of the data buffer that is filled.
     */
    double_t data_load_factor() const {
        return static_cast<double_t>(data_size_) / max_data_size_;
    }

    /**
     * @brief The percentage of the offset and validity buffers that is filled if any.
     */
    double_t cell_load_factor() const {
        return static_cast<double_t>(num_cells_) / max_num_cells_;
    }

   protected:
    size_t num_cells_;

    // Data size which is calculated different for var vs non-var
    size_t data_size_;

    size_t max_num_cells_;

    // Data size which is calculated different for var vs non-var
    size_t max_data_size_;

   private:
    void attach_buffer(Query& query);

    void attach_subarray(Subarray& subarray) const;

    template <typename T>
    void attach_range(Subarray& subarray) const {
        subarray.add_range(name_, data<T>()[0], data<T>()[1]);
    }

    //===================================================================
    //= private non-static
    //===================================================================

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
    std::optional<Enumeration> enumeration_;

    // True if the array has at least one enumerations
    bool has_enumeration_ = false;

    // Enumerations (optional)
    std::vector<std::string> enums_;

    // Enumerations (optional) as string and offsets
    std::string enum_str_;
    std::vector<uint32_t> enum_offsets_;
    bool is_ordered_ = false;
};

class ReadColumnBuffer : public ColumnBuffer {
   public:
    using ColumnBuffer::ColumnBuffer;

    virtual ~ReadColumnBuffer();

    /**
     * @brief Size num_cells_ to match the read query results.
     *
     * @param query TileDB query
     */
    size_t update_size(const Query& query);

    /**
     * @brief Return a view of the ColumnBuffer data.
     *
     * @tparam T Data type
     * @return std::span<T> data view
     */
    template <typename T>
    std::span<T> data() {
        return std::span<T>(reinterpret_cast<T*>(data().data()), size());
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

    /**
     * @brief Convert the data bytemap to a bitmap in place.
     *
     */
    void data_to_bitmap() {
        ColumnBuffer::to_bitmap(data<uint8_t>());
    }

    /**
     * @brief Convert the validity bytemap to a bitmap in place.
     *
     */
    void validity_to_bitmap() {
        ColumnBuffer::to_bitmap(validity());
    }

    void validity_to_bitmap(std::span<uint8_t> bitmap) const {
        ColumnBuffer::to_bitmap(validity(), bitmap);
    }

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
        size_t num_cells,
        size_t num_bytes,
        bool is_var = false,
        bool is_nullable = false,
        std::optional<Enumeration> enumeration = std::nullopt,
        bool is_ordered = false);

    CArrayColumnBuffer() = delete;
    CArrayColumnBuffer(const CArrayColumnBuffer&) = delete;
    CArrayColumnBuffer(CArrayColumnBuffer&&) = default;

    ~CArrayColumnBuffer();

    std::span<std::byte> data() override;
    std::span<uint64_t> offsets() override;
    std::span<uint8_t> validity() override;

    std::span<const std::byte> data() const override;
    std::span<const uint64_t> offsets() const override;
    std::span<const uint8_t> validity() const override;

   private:
    std::unique_ptr<std::byte[]> data_;
    std::unique_ptr<uint64_t[]> offsets_;
    std::unique_ptr<uint8_t[]> validity_;
};

/**
 * @brief Class to store data for a TileDB dimension or attribute.
 *
 */
class VectorColumnBuffer : public ReadColumnBuffer {
    // This "medium size" -- 1 GiB -- is a good balance between improved remote
    // I/O performance (bigger = better) and friendliness for smaller hardware
    // (e.g.  tiny CI runners). This value is good for general in-between
    // hardware including modern laptops and a broad range of EC2 instances. CI
    // can ask for smaller; power-server users can ask for larger.
    inline static const size_t DEFAULT_ALLOC_BYTES = 1 << 30;
    inline static const std::string CONFIG_KEY_INIT_BYTES = "soma.init_buffer_bytes";

   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a ColumnBuffer from an array and column name.
     *
     * @param array TileDB array
     * @param name TileDB dimension or attribute name
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> create(std::shared_ptr<Array> array, std::string_view name);

    //===================================================================
    //= public non-static
    //===================================================================

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
    VectorColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        size_t num_bytes,
        bool is_var = false,
        bool is_nullable = false,
        std::optional<Enumeration> enumeration = std::nullopt,
        bool is_ordered = false);

    VectorColumnBuffer() = delete;
    VectorColumnBuffer(const VectorColumnBuffer&) = delete;
    VectorColumnBuffer(VectorColumnBuffer&&) = default;

    ~VectorColumnBuffer();

    std::span<const std::byte> data() const override;
    std::span<const uint64_t> offsets() const override;
    std::span<const uint8_t> validity() const override;

    std::span<std::byte> data() override;
    std::span<uint64_t> offsets() override;
    std::span<uint8_t> validity() override;

   private:
    //===================================================================
    //= private static
    //===================================================================

    /**
     * @brief Allocate and return a ColumnBuffer.
     *
     * @param config TileDB Config
     * @param name Column name
     * @param type TileDB datatype
     * @param is_var True if variable length data
     * @param is_nullable True if nullable data
     * @param enumeration Optional Enumeration associated with column
     * @param is_ordered Optional Enumeration is ordered
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> alloc(
        Config config,
        std::string_view name,
        tiledb_datatype_t type,
        bool is_var,
        bool is_nullable,
        std::optional<Enumeration> enumeration,
        bool is_ordered);

    //===================================================================
    //= private non-static
    //===================================================================

    // Data buffer.
    std::vector<std::byte> data_;

    // Offsets buffer (optional).
    std::vector<uint64_t> offsets_;

    // Validity buffer (optional).
    std::vector<uint8_t> validity_;
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
        std::shared_ptr<Array> array,
        std::string_view name,
        size_t num_elements,
        DataStorage data_buffer,
        OffsetStorage offsets_buffer,
        std::unique_ptr<uint8_t[]> validity_buffer,
        bool copy_buffers = false) {
        if (array->query_type() != TILEDB_WRITE) {
            throw TileDBSOMAError("[WriteColumnBuffer] WriteColumnBuffer should only be used during write.");
        }

        auto schema = array->schema();
        auto name_str = std::string(name);  // string for TileDB API

        if (schema.has_attribute(name_str)) {
            auto attr = schema.attribute(name_str);
            auto type = attr.type();
            bool is_var = attr.cell_val_num() == TILEDB_VAR_NUM;
            bool is_nullable = attr.nullable();

            if (!is_var && attr.cell_val_num() != 1) {
                throw TileDBSOMAError("[WriteColumnBuffer] Values per cell > 1 is not supported: " + name_str);
            }

            size_t num_bytes = is_var ? offsets_buffer[num_elements] : num_elements * impl::type_size(type);

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
                throw TileDBSOMAError("[WriteColumnBuffer] Values per cell > 1 is not supported: " + name_str);
            }

            size_t num_bytes = is_var ? offsets_buffer[num_elements] : num_elements * impl::type_size(type);

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

        throw TileDBSOMAError("[WriteColumnBuffer] Column name not found: " + name_str);
    }

    template <typename DataStorage, typename OffsetStorage>
        requires is_data_buffer<DataStorage> && is_offset_buffer<OffsetStorage>
    WriteColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        size_t num_bytes,
        bool is_var,
        bool is_nullable,
        DataStorage data_buffer,
        OffsetStorage offsets_buffer,
        std::unique_ptr<uint8_t[]> validity_buffer,
        bool copy_buffers = false)
        : ColumnBuffer(
              name, type, num_cells, num_cells, num_bytes, num_bytes, is_var, is_nullable, std::nullopt, false) {
        if constexpr (std::is_same_v<std::unique_ptr<std::byte[]>, DataStorage>) {
            // Data buffer ownership is passed to the write column buffer
            data_buffer_ = std::move(data_buffer);
            data_ = std::span<const std::byte>(data_buffer_.get(), data_size_);
        } else {
            if (copy_buffers) {
                data_buffer_ = std::make_unique_for_overwrite<std::byte[]>(data_size_);
                std::memcpy(data_buffer_.get(), data_buffer, data_size_);

                data_ = std::span<const std::byte>(data_buffer_.get(), data_size_);
            } else {
                data_ = std::span<const std::byte>(reinterpret_cast<const std::byte*>(data_buffer), data_size_);
            }
        }

        if (is_var) {
            if constexpr (std::is_same_v<std::unique_ptr<uint64_t[]>, OffsetStorage>) {
                // Offset buffer ownership is passed to the write column buffer
                offsets_buffer_ = std::move(offsets_buffer);
                offsets_ = std::span<const uint64_t>(offsets_buffer_.get(), num_cells_ + 1);
            } else if constexpr (std::is_same_v<uint32_t, std::remove_const_t<std::remove_pointer_t<OffsetStorage>>>) {
                // A new offset buffer will be contructed using the input offset buffer
                offsets_buffer_ = std::make_unique_for_overwrite<uint64_t[]>(num_cells_ + 1);

                // Fill the buffer with offset and cast to the correct type
                std::transform(
                    offsets_buffer, offsets_buffer + num_cells_ + 1, offsets_buffer_.get(), [](const uint32_t& offset) {
                        return static_cast<uint64_t>(offset);
                    });

                offsets_ = std::span<const uint64_t>(offsets_buffer_.get(), num_cells_ + 1);
            } else {
                if (copy_buffers) {
                    offsets_buffer_ = std::make_unique_for_overwrite<uint64_t[]>(num_cells_ + 1);
                    std::memcpy(offsets_buffer_.get(), offsets_buffer, (num_cells_ + 1) * sizeof(uint64_t));

                    offsets_ = std::span<const uint64_t>(offsets_buffer_.get(), num_cells_ + 1);
                } else {
                    offsets_ = std::span<const uint64_t>(offsets_buffer, num_cells_ + 1);
                }
            }
        }

        if (is_nullable) {
            if (validity_buffer) {
                validity_buffer_ = std::move(validity_buffer);
                validity_ = std::span<const uint8_t>(validity_buffer_.get(), num_cells);
            } else {
                validity_buffer_ = std::make_unique_for_overwrite<uint8_t[]>(num_cells);
                std::fill(validity_buffer_.get(), validity_buffer_.get() + num_cells, 1);
                validity_ = std::span<const uint8_t>(validity_buffer_.get(), num_cells);
            }
        } else {
            // As of version 1.15.6 we were throwing here. However, we found a
            // compatibility issue with pyarrow versions below 17. Thus we log and
            // continue.
            if (validity_buffer) {
                LOG_DEBUG(
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
    std::span<const std::byte> data_;
    std::span<const uint64_t> offsets_;
    std::span<const uint8_t> validity_;

    std::unique_ptr<std::byte[]> data_buffer_;
    std::unique_ptr<uint64_t[]> offsets_buffer_;
    std::unique_ptr<uint8_t[]> validity_buffer_;
};

}  // namespace tiledbsoma
#endif
