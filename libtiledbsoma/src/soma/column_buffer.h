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

#include <span>
#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include "../utils/arrow_adapter.h"
#include "../utils/common.h"
#include "soma_context.h"

namespace tiledbsoma {

using namespace tiledb;

/**
 * @brief Class to store data for a TileDB dimension or attribute.
 *
 */
class ColumnBuffer {
    // This "medium size" -- 1 GiB -- is a good balance between improved remote
    // I/O performance (bigger = better) and friendliness for smaller hardware
    // (e.g.  tiny CI runners). This value is good for general in-between
    // hardware including modern laptops and a broad range of EC2 instances. CI
    // can ask for smaller; power-server users can ask for larger.
    inline static const size_t DEFAULT_ALLOC_BYTES = 1 << 30;
    inline static const std::string
        CONFIG_KEY_INIT_BYTES = "soma.init_buffer_bytes";

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
    static std::shared_ptr<ColumnBuffer> create(
        std::shared_ptr<Array> array, std::string_view name);

    /**
     * @brief Convert a bytemap to a bitmap in place.
     *
     */
    static void to_bitmap(std::span<uint8_t> bytemap);

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
    ColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        size_t num_bytes,
        bool is_var = false,
        bool is_nullable = false,
        std::optional<Enumeration> enumeration = std::nullopt,
        bool is_ordered = false);

    ColumnBuffer() = delete;
    ColumnBuffer(const ColumnBuffer&) = delete;
    ColumnBuffer(ColumnBuffer&&) = default;

    ~ColumnBuffer();

    /**
     * @brief Attach this ColumnBuffer to a TileDB query.
     *
     * @param query TileDB query
     */
    void attach(Query& query, std::optional<Subarray> subarray = std::nullopt);

    template <typename T>
    void set_data(
        uint64_t num_elems,
        const void* data,
        T* offsets,
        const std::optional<std::vector<uint8_t>>& validity = std::nullopt) {
        num_cells_ = num_elems;

        // Ensure the offset type is either uint32_t* or uint64_t*
        static_assert(
            std::is_same_v<T, uint32_t> || std::is_same_v<T, uint64_t>,
            "offsets must be either uint32_t* or uint64_t*");

        if (offsets != nullptr) {
            offsets_ = std::vector<uint64_t>(
                (T*)offsets, (T*)offsets + num_elems + 1);

            data_size_ = offsets_[num_elems];
            data_.assign((std::byte*)data, (std::byte*)data + data_size_);
        } else {
            data_size_ = num_elems;
            data_.assign(
                (std::byte*)data, (std::byte*)data + num_elems * type_size_);
        }

        if (validity.has_value()) {
            validity_ = *validity;
        } else {
            if (is_nullable_) {
                validity_.assign(num_elems, 1);  // Default all to valid (1)
            }
        }
    }

    /**
     * @brief Size num_cells_ to match the read query results.
     *
     * @param query TileDB query
     */
    size_t update_size(const Query& query);

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
    uint64_t data_size() {
        return data_size_;
    }

    /**
     * @brief Return a view of the ColumnBuffer data.
     *
     * @tparam T Data type
     * @return std::span<T> data view
     */
    template <typename T>
    std::span<T> data() {
        return std::span<T>((T*)data_.data(), num_cells_);
    }

    /**
     * @brief Return data in a vector of binary buffers.
     *
     * @return std::vector<std::vector<std::byte>>
     */
    std::vector<std::vector<std::byte>> binaries();

    /**
     * @brief Return data in a vector of strings.
     *
     * @return std::vector<std::string>
     */
    std::vector<std::string> strings();

    /**
     * @brief Return a string_view of the string at the provided cell index.
     *
     * @param index Cell index
     * @return std::string_view string view
     */
    std::string_view string_view(uint64_t index);

    /**
     * @brief Return a view of the ColumnBuffer offsets.
     *
     * @return std::span<uint64_t> offsets view
     */
    std::span<uint64_t> offsets() {
        if (!is_var_) {
            throw TileDBSOMAError(
                "[ColumnBuffer] Offsets buffer not defined for " + name_);
        }

        return std::span<uint64_t>(offsets_.data(), num_cells_);
    }

    /**
     * @brief Return a view of the validity buffer.
     *
     * @return std::span<uint8_t> validity view
     */
    std::span<uint8_t> validity() {
        if (!is_nullable_) {
            throw TileDBSOMAError(
                "[ColumnBuffer] Validity buffer not defined for " + name_);
        }
        return std::span<uint8_t>(validity_.data(), num_cells_);
    }

    /**
     * @brief Return the name of the buffer.
     *
     * @return std::string_view
     */
    std::string_view name() {
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
            throw TileDBSOMAError(
                "[ColumnBuffer] No enumeration defined for " + name_);
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
            throw TileDBSOMAError(
                "[ColumnBuffer] No enumeration defined for " + name_);
        }
        return std::span<uint32_t>(enum_offsets_.data(), enum_offsets_.size());
    }

    /**
     * @brief Return optional enumeration string
     *
     */
    std::span<char> enum_string() {
        if (!has_enumeration_) {
            throw TileDBSOMAError(
                "[ColumnBuffer] No enumeration defined for " + name_);
        }
        return std::span<char>(enum_str_.data(), enum_str_.length());
    }

    /**
     * @brief Return true if the buffer contains an ordered enumeration.
     */
    bool is_ordered() const {
        return is_ordered_;
    }

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

    void attach_buffer(Query& query);

    void attach_subarray(Subarray& subarray);

    template <typename T>
    void attach_range(Subarray& subarray) {
        subarray.add_range(name_, data<T>()[0], data<T>()[1]);
    }

    //===================================================================
    //= private non-static
    //===================================================================

    // Name of the column from the schema.
    std::string name_;

    // Data type of the column from the schema.
    tiledb_datatype_t type_;

    // Data size which is calculated different for var vs non-var
    uint64_t data_size_;

    // Bytes per element.
    uint64_t type_size_;

    // Number of cells.
    uint64_t num_cells_;

    // If true, the data type is variable length
    bool is_var_;

    // If true, the data is nullable
    bool is_nullable_;

    // If applicable, the Enumeration associated with the column
    std::optional<Enumeration> enumeration_;

    // Data buffer.
    std::vector<std::byte> data_;

    // Offsets buffer (optional).
    std::vector<uint64_t> offsets_;

    // Validity buffer (optional).
    std::vector<uint8_t> validity_;

    // True if the array has at least one enumerations
    bool has_enumeration_ = false;

    // Enumerations (optional)
    std::vector<std::string> enums_;

    // Enumerations (optional) as string and offsets
    std::string enum_str_;
    std::vector<uint32_t> enum_offsets_;
    bool is_ordered_ = false;
};

}  // namespace tiledbsoma
#endif
