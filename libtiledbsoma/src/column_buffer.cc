/**
 * @file   column_buffer.cc
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
 * This file defines the a ColumBuffer class.
 */

#include "tiledbsoma/column_buffer.h"
#include "tiledbsoma/array_buffers.h"
#include "tiledbsoma/common.h"
#include "tiledbsoma/logger_public.h"

namespace tiledbsoma {

using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::shared_ptr<ColumnBuffer> ColumnBuffer::create(
    std::shared_ptr<Array> array, std::string_view name) {
    auto name_str = std::string(name);  // string for TileDB API
    auto schema = array->schema();

    if (schema.has_attribute(name_str)) {
        auto attr = schema.attribute(name_str);
        bool is_var = attr.cell_val_num() == TILEDB_VAR_NUM;
        bool is_nullable = attr.nullable();

        if (!is_var && attr.cell_val_num() != 1) {
            throw TileDBSOMAError(
                "[ColumnBuffer] Values per cell > 1 is not supported: " +
                name_str);
        }

        return ColumnBuffer::alloc(
            array, attr.name(), attr.type(), is_var, is_nullable);

    } else if (schema.domain().has_dimension(name_str)) {
        auto dim = schema.domain().dimension(name_str);
        bool is_var = dim.cell_val_num() == TILEDB_VAR_NUM ||
                      dim.type() == TILEDB_STRING_ASCII ||
                      dim.type() == TILEDB_STRING_UTF8;

        if (!is_var && dim.cell_val_num() != 1) {
            throw TileDBSOMAError(
                "[ColumnBuffer] Values per cell > 1 is not supported: " +
                name_str);
        }

        return ColumnBuffer::alloc(
            array, dim.name(), dim.type(), is_var, false);
    }

    throw TileDBSOMAError("[ColumnBuffer] Column name not found: " + name_str);
}

void ColumnBuffer::to_bitmap(tcb::span<uint8_t> bytemap) {
    int i_dst = 0;
    for (unsigned int i_src = 0; i_src < bytemap.size(); i_src++) {
        // Overwrite every 8 bytes with a one-byte bitmap
        if (i_src % 8 == 0) {
            // Each bit in the bitmap corresponds to one byte in the bytemap
            // Note: the bitmap must be byte-aligned (8 bits)
            int bitmap = 0;
            for (unsigned int i = i_src; i < i_src + 8 && i < bytemap.size();
                 i++) {
                bitmap |= bytemap[i] << (i % 8);
            }
            bytemap[i_dst++] = bitmap;
        }
    }
}

//===================================================================
//= public non-static
//===================================================================

ColumnBuffer::ColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    size_t num_cells,
    size_t num_bytes,
    bool is_var,
    bool is_nullable)
    : name_(name)
    , type_(type)
    , type_size_(tiledb::impl::type_size(type))
    , num_cells_(0)
    , is_var_(is_var)
    , is_nullable_(is_nullable) {
    LOG_DEBUG(fmt::format(
        "[ColumnBuffer] '{}' {} bytes is_var={} is_nullable={}",
        name,
        num_bytes,
        is_var_,
        is_nullable_));
    // Call reserve() to allocate memory without initializing the contents.
    // This reduce the time to allocate the buffer and reduces the
    // resident memory footprint of the buffer.
    data_.reserve(num_bytes);
    if (is_var_) {
        offsets_.reserve(num_cells + 1);  // extra offset for arrow
    }
    if (is_nullable_) {
        validity_.reserve(num_cells);
    }
}

void ColumnBuffer::attach(Query& query) {
    // We cannot use:
    // `set_data_buffer(const std::string& name, std::vector<T>& buf)`
    // because data_ is allocated with reserve() and data_.size()
    // does not represent the actual size of the buffer.
    query.set_data_buffer(
        name_, (void*)data_.data(), data_.capacity() / type_size_);
    if (is_var_) {
        // Remove one offset for TileDB, which checks that the
        // offsets and validity buffers are the same size
        query.set_offsets_buffer(
            name_, offsets_.data(), offsets_.capacity() - 1);
    }
    if (is_nullable_) {
        query.set_validity_buffer(
            name_, validity_.data(), validity_.capacity());
    }
}

size_t ColumnBuffer::update_size(const Query& query) {
    auto [num_offsets, num_elements] = query.result_buffer_elements()[name_];

    if (is_var()) {
        num_cells_ = num_offsets;
        // Set the extra offset value for arrow.
        offsets_[num_offsets] = num_elements;
    } else {
        num_cells_ = num_elements;
    }

    return num_cells_;
}

std::vector<std::string> ColumnBuffer::strings() {
    std::vector<std::string> result;

    for (size_t i = 0; i < num_cells_; i++) {
        result.emplace_back(std::string(string_view(i)));
    }

    return result;
}

std::string_view ColumnBuffer::string_view(uint64_t index) {
    auto start = offsets_[index];
    auto len = offsets_[index + 1] - start;
    return std::string_view((char*)(data_.data() + start), len);
}

//===================================================================
//= private static
//===================================================================

std::shared_ptr<ColumnBuffer> ColumnBuffer::alloc(
    std::shared_ptr<Array> array,
    std::string_view name,
    tiledb_datatype_t type,
    bool is_var,
    bool is_nullable) {
    // Set number of bytes for the data buffer. Override with a value from
    // the config if present.
    auto num_bytes = DEFAULT_ALLOC_BYTES;
    auto config = array->schema().context().config();
    if (config.contains(CONFIG_KEY_INIT_BYTES)) {
        auto value_str = config.get(CONFIG_KEY_INIT_BYTES);
        try {
            num_bytes = std::stoull(value_str);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(fmt::format(
                "[ColumnBuffer] Error parsing {}: '{}' ({})",
                CONFIG_KEY_INIT_BYTES,
                value_str,
                e.what()));
        }
    }

    bool is_dense = array->schema().array_type() == TILEDB_DENSE;
    if (is_dense) {
        // TODO: Handle dense arrays similar to tiledb python module
    }

    // For variable length column types, allocate an extra num_bytes to hold
    //   offset values. The number of cells is the set by the size of the
    //   offset type.
    // For non-variable length column types, the number of cells is computed
    //   from the type size.
    size_t num_cells = is_var ? num_bytes / sizeof(uint64_t) :
                                num_bytes / tiledb::impl::type_size(type);

    return std::make_shared<ColumnBuffer>(
        name, type, num_cells, num_bytes, is_var, is_nullable);
}

}  // namespace tiledbsoma
