/**
 * @file   column_buffer.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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

#include <span/span.hpp>
#include <tiledb/tiledb>

#include "tiledbsoma/common.h"
#include "tiledbsoma/logger_public.h"

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
     */
    ColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        size_t num_bytes,
        bool is_var = false,
        bool is_nullable = false);

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
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> alloc(
        std::shared_ptr<Array> array,
        std::string_view name,
        tiledb_datatype_t type,
        bool is_var,
        bool is_nullable);

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

    // Data buffer.
    std::vector<std::byte> data_;

    // Offsets buffer (optional).
    std::vector<uint64_t> offsets_;

    // Validity buffer (optional).
    std::vector<uint8_t> validity_;
};

}  // namespace tiledbsoma
#endif
