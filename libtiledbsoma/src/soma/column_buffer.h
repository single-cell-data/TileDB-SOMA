/**
 * @file   column_buffer.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022-2023 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 *   This declares the column buffer API
 */

#ifndef COLUMN_BUFFER_H
#define COLUMN_BUFFER_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include "../utils/common.h"
#include "../utils/logger.h"
#include "span/span.hpp"

namespace tiledbsoma {

using namespace tiledb;

/**
 * @brief Class to store data for a TileDB dimension or attribute.
 *
 */
class ColumnBuffer {
    inline static const size_t DEFAULT_ALLOC_BYTES = 1 << 24;  // 16 MiB
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
     * @brief Create a ColumnBuffer from a schema, column name, and data.
     *
     * @param array TileDB array
     * @param name TileDB dimension or attribute name
     * @param data Data to set in buffer
     * @return ColumnBuffer
     */
    template <typename T>
    static std::shared_ptr<ColumnBuffer> create(
        std::shared_ptr<Array> array,
        std::string_view name,
        std::vector<T> data) {
        auto column_buff = ColumnBuffer::create(array, name);
        column_buff->num_cells_ = data.size();
        column_buff->data_.resize(data.size());
        column_buff->data_.assign(
            reinterpret_cast<std::byte*>(data.data()),
            reinterpret_cast<std::byte*>(data.data() + data.size()));
        return column_buff;
    }

    /**
     * @brief Convert a bytemap to a bitmap in place.
     *
     */
    static void to_bitmap(tcb::span<uint8_t> bytemap);

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

    ~ColumnBuffer() {
        LOG_TRACE(fmt::format("[ColumnBuffer] release '{}'", name_));
    }

    /**
     * @brief Attach this ColumnBuffer to a TileDB query.
     *
     * @param query TileDB query
     */
    void attach(Query& query);

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
     * @brief Return a view of the ColumnBuffer data.
     *
     * @tparam T Data type
     * @return tcb::span<T> data view
     */
    template <typename T>
    tcb::span<T> data() {
        return tcb::span<T>((T*)data_.data(), num_cells_);
    }

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
     * @return tcb::span<uint64_t> offsets view
     */
    tcb::span<uint64_t> offsets() {
        if (!is_var_) {
            throw TileDBSOMAError(
                "[ColumnBuffer] Offsets buffer not defined for " + name_);
        }

        return tcb::span<uint64_t>(offsets_.data(), num_cells_);
    }

    /**
     * @brief Return a view of the validity buffer.
     *
     * @return tcb::span<uint8_t> validity view
     */
    tcb::span<uint8_t> validity() {
        if (!is_nullable_) {
            throw TileDBSOMAError(
                "[ColumnBuffer] Validity buffer not defined for " + name_);
        }
        return tcb::span<uint8_t>(validity_.data(), num_cells_);
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

    std::optional<Enumeration> get_enumeration() const {
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
            cumlen += s.length();
        }
        enum_offsets_[n_vec] = cumlen;
    }

    /**
     * @brief Return optional enumeration offsets vector
     *
     */
    tcb::span<uint32_t> enum_offsets() {
        if (!has_enumeration_) {
            throw TileDBSOMAError(
                "[ColumnBuffer] No enumeration defined for " + name_);
        }
        return tcb::span<uint32_t>(enum_offsets_.data(), enum_offsets_.size());
    }

    /**
     * @brief Return optional enumeration string
     *
     */
    tcb::span<char> enum_string() {
        if (!has_enumeration_) {
            throw TileDBSOMAError(
                "[ColumnBuffer] No enumeration defined for " + name_);
        }
        return tcb::span<char>(enum_str_.data(), enum_str_.length());
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
     * @param array TileDB array
     * @param name Column name
     * @param type TileDB datatype
     * @param is_var True if variable length data
     * @param is_nullable True if nullable data
     * @param enumeration Optional Enumeration associated with column
     * @param is_ordered Optional Enumeration is ordered
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> alloc(
        ArraySchema schema,
        std::string_view name,
        tiledb_datatype_t type,
        bool is_var,
        bool is_nullable,
        std::optional<Enumeration> enumeration,
        bool is_ordered);

    //===================================================================
    //= private non-static
    //===================================================================

    // Name of the column from the schema.
    std::string name_;

    // Data type of the column from the schema.
    tiledb_datatype_t type_;

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
