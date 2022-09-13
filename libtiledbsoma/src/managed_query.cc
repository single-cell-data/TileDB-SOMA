/**
 * @file   managed_query.cc
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
 * This file defines the performing TileDB queries.
 */

#include "tiledbsoma/managed_query.h"
#include "tiledbsoma/common.h"
#include "tiledbsoma/logger_public.h"

namespace tiledbsoma {

using namespace tiledb;

//===================================================================
//= public non-static
//===================================================================

ManagedQuery::ManagedQuery(std::shared_ptr<Array> array, std::string name)
    : array_(array)
    , name_(name)
    , schema_(array->schema()) {
    query_ = std::make_unique<Query>(schema_.context(), *array);
    subarray_ = std::make_unique<Subarray>(schema_.context(), *array);

    if (array->schema().array_type() == TILEDB_SPARSE) {
        query_->set_layout(TILEDB_UNORDERED);
    } else {
        query_->set_layout(TILEDB_ROW_MAJOR);
    }
}

void ManagedQuery::select_columns(
    std::vector<std::string> names, bool if_not_empty) {
    // Return if we are selecting all columns (columns_ is empty) and we want to
    // continue selecting all columns (if_not_empty == true).
    if (if_not_empty && columns_.empty()) {
        return;
    }

    for (auto& name : names) {
        // Name is not an attribute or dimension.
        if (!schema_.has_attribute(name) &&
            !schema_.domain().has_dimension(name)) {
            LOG_DEBUG(fmt::format(
                "[ManagedQuery] [{}] Invalid column selected: {}",
                name_,
                name));
            invalid_columns_selected_ = true;
        } else {
            columns_.insert(name);
        }
    }
}

size_t ManagedQuery::submit() {
    auto status = query_->query_status();

    // Query is complete, return 0 cells read
    if (status == Query::Status::COMPLETE) {
        return 0;
    }

    // Query is uninitialized, allocate and attach buffers
    if (status == Query::Status::UNINITIALIZED) {
        // Set the subarray for range slicing
        query_->set_subarray(*subarray_);

        // If no columns were selected, select all columns
        if (!columns_.size()) {
            for (const auto& dim : array_->schema().domain().dimensions()) {
                columns_.insert(dim.name());
            }
            for (const auto& [name, attr] : array_->schema().attributes()) {
                (void)attr;
                columns_.insert(name);
            }
        }

        // Allocate and attach buffers
        for (auto& name : columns_) {
            LOG_DEBUG(fmt::format(
                "[ManagedQuery] [{}] Adding buffer for column '{}'",
                name_,
                name));
            buffers_.emplace(name, ColumnBuffer::create(array_, name));
            buffers_[name]->attach(*query_);
        }
    }

    // Submit query
    LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Submit query", name_));
    query_->submit();
    status = query_->query_status();
    LOG_DEBUG(fmt::format(
        "[ManagedQuery] [{}] Query status = {}", name_, (int)status));

    // If the query was ever incomplete, the result buffers contents are not
    // complete.
    if (status == Query::Status::INCOMPLETE) {
        results_complete_ = false;
    }

    // Update ColumnBuffer size to match query results
    size_t num_cells = 0;
    for (auto& [name, buffer] : buffers_) {
        num_cells = buffer->update_size(*query_);
        LOG_DEBUG(fmt::format(
            "[ManagedQuery] [{}] Buffer {} cells={}", name_, name, num_cells));
    }
    total_num_cells_ += num_cells;

    // TODO: retry the query with larger buffers
    if (status == Query::Status::INCOMPLETE && !num_cells) {
        throw TileDBSCError(
            fmt::format("[ManagedQuery] [{}] Buffers are too small.", name_));
    }

    return num_cells;
}

};  // namespace tiledbsoma
