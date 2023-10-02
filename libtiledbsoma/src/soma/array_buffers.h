/**
 * @file   array_buffers.h
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
 *   This declares the array buffers API
 */

#ifndef ARRAY_BUFFERS_H
#define ARRAY_BUFFERS_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <tiledb/tiledb>

#include "../utils/arrow_adapter.h"
#include "../utils/common.h"
#include "../utils/logger.h"
#include "column_buffer.h"

namespace tiledbsoma {

using namespace tiledb;

class ArrayBuffers {
   public:
    ArrayBuffers() = default;
    ArrayBuffers(const ArrayBuffers&) = delete;
    ArrayBuffers(ArrayBuffers&&) = default;
    ~ArrayBuffers() = default;

    static std::shared_ptr<ArrayBuffers> from_arrow(
        ArraySchema& schema,
        ArrowSchema& arrow_schema,
        ArrowArray& arrow_array) {
        auto array_buffer = std::make_shared<ArrayBuffers>();
        for (int64_t i = 0; i < arrow_array.n_children; ++i) {
            auto name = arrow_schema.children[i]->name;
            auto data = arrow_array.children[i]->buffers[1];
            auto typeinfo = ArrowAdapter::arrow_type_to_tiledb(
                arrow_schema.children[i]);
            array_buffer->emplace(
                name,
                ColumnBuffer::create(
                    schema,
                    name,
                    data,
                    arrow_array.children[i]->length,
                    typeinfo.elem_size));
        }
        return array_buffer;
    }

    /**
     * @brief Return the buffer with the given name.
     *
     * @param name Column name
     * @return std::shared_ptr<ColumnBuffer> Column buffer
     */
    std::shared_ptr<ColumnBuffer> at(const std::string& name) {
        if (!contains(name)) {
            throw TileDBSOMAError(
                fmt::format("[ArrayBuffers] column '{}' does not exist", name));
        }
        return buffers_[name];
    }

    /**
     * @brief Return true if a buffer with the given name exists.
     *
     * @param name Column name
     * @return True if a buffer with the given name exists
     */
    bool contains(const std::string& name) {
        return buffers_.find(name) != buffers_.end();
    }

    /**
     * @brief Add a column buffer with the given name to the ArrayBuffers,
     * maintaining the insertion order.
     *
     * @param name Column name
     * @param buffer Column buffer
     */
    void emplace(
        const std::string& name, std::shared_ptr<ColumnBuffer> buffer) {
        if (contains(name)) {
            throw TileDBSOMAError(
                fmt::format("[ArrayBuffers] column '{}' already exists", name));
        }
        names_.push_back(name);
        buffers_.emplace(name, buffer);
    }

    /**
     * @brief Returns the ordered vector of names.
     *
     * @return const std::vector<std::string>& Vector of names
     */
    const std::vector<std::string>& names() {
        return names_;
    }

    /**
     * @brief Returns the number of rows in the array buffer.
     *
     * @return uint64_t Number of rows
     */
    uint64_t num_rows() const {
        return buffers_.at(names_.front())->size();
    }

   private:
    // A vector of column names that maintains the order the columns were added
    std::vector<std::string> names_;

    // Map: column name -> ColumnBuffer
    std::unordered_map<std::string, std::shared_ptr<ColumnBuffer>> buffers_;
};

}  // namespace tiledbsoma

#endif
