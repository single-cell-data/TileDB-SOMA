/**
 * @file   array_buffers.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This declares the array buffers API
 */

#ifndef ARRAY_BUFFERS_H
#define ARRAY_BUFFERS_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <tiledb/tiledb>

#include "../utils/common.h"
#include "column_buffer.h"

namespace tiledbsoma {

using namespace tiledb;

class ArrayBuffers {
   public:
    ArrayBuffers() = default;
    ArrayBuffers(const ArrayBuffers&) = default;
    ArrayBuffers(ArrayBuffers&&) = default;
    ~ArrayBuffers() = default;

    /**
     * @brief Return the buffer with the given name.
     *
     * @param name Column name
     * @return std::shared_ptr<ColumnBuffer> Column buffer
     */
    std::shared_ptr<ColumnBuffer> at(const std::string& name);

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
    void emplace(const std::string& name, std::shared_ptr<ColumnBuffer> buffer);

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
