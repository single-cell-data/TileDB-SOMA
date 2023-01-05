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

ManagedQuery::ManagedQuery(std::shared_ptr<Array> array, std::string_view name)
    : array_(array)
    , name_(name)
    , schema_(std::make_shared<ArraySchema>(array->schema())) {
    reset();
}

void ManagedQuery::reset() {
    query_ = std::make_unique<Query>(schema_->context(), *array_);
    subarray_ = std::make_unique<Subarray>(schema_->context(), *array_);

    if (array_->schema().array_type() == TILEDB_SPARSE) {
        query_->set_layout(TILEDB_UNORDERED);
    } else {
        query_->set_layout(TILEDB_ROW_MAJOR);
    }

    subarray_range_set_ = false;
    subarray_range_empty_ = true;
    columns_.clear();
    results_complete_ = true;
    total_num_cells_ = 0;
    buffers_.reset();
    query_submitted_ = false;
}

void ManagedQuery::select_columns(
    const std::vector<std::string>& names, bool if_not_empty) {
    // Return if we are selecting all columns (columns_ is empty) and we want to
    // continue selecting all columns (if_not_empty == true).
    if (if_not_empty && columns_.empty()) {
        return;
    }

    for (auto& name : names) {
        // Name is not an attribute or dimension.
        if (!schema_->has_attribute(name) &&
            !schema_->domain().has_dimension(name)) {
            LOG_WARN(fmt::format(
                "[TileDB-SOMA::ManagedQuery] [{}] Invalid column selected: {}",
                name_,
                name));
        } else {
            columns_.push_back(name);
        }
    }
}

void ManagedQuery::submit() {
    // Throw error if submit is called again before reading the results
    if (query_submitted_) {
        throw TileDBSOMAError(fmt::format(
            "[ManagedQuery][{}] read results before calling submit again",
            name_));
    }

    // If the query is complete, return so we do not submit it again
    auto status = query_->query_status();
    if (status == Query::Status::COMPLETE) {
        return;
    }

    // If the query is uninitialized, set the subarray for the query
    if (status == Query::Status::UNINITIALIZED) {
        // Dense array must have a subarray set. If the array is dense and no
        // ranges have been set, add a range for the array's entire non-empty
        // domain on dimension 0.
        if (array_->schema().array_type() == TILEDB_DENSE &&
            !subarray_range_set_) {
            auto non_empty_domain = array_->non_empty_domain<int64_t>(0);
            subarray_->add_range(
                0, non_empty_domain.first, non_empty_domain.second);

            LOG_DEBUG(fmt::format(
                "[ManagedQuery] Add full NED range to dense subarray = (0, {}, "
                "{})",
                non_empty_domain.first,
                non_empty_domain.second));
        }

        // Set the subarray for range slicing
        query_->set_subarray(*subarray_);
    }

    // If no columns were selected, select all columns.
    // Add dims and attrs in the same order as specified in the schema
    if (columns_.empty()) {
        for (const auto& dim : array_->schema().domain().dimensions()) {
            columns_.push_back(dim.name());
        }
        int attribute_num = array_->schema().attribute_num();
        for (int i = 0; i < attribute_num; i++) {
            columns_.push_back(array_->schema().attribute(i).name());
        }
    }

    // Allocate and attach buffers
    LOG_TRACE("[ManagedQuery] allocate new buffers");
    buffers_ = std::make_shared<ArrayBuffers>();
    for (auto& name : columns_) {
        LOG_DEBUG(fmt::format(
            "[ManagedQuery] [{}] Adding buffer for column '{}'", name_, name));
        buffers_->emplace(name, ColumnBuffer::create(array_, name));
        buffers_->at(name)->attach(*query_);
    }

    // Submit query
    LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Submit query", name_));

    // Do not submit if the query contains only empty ranges
    if (!is_empty_query()) {
        query_->submit();
    }
    query_submitted_ = true;
}

std::shared_ptr<ArrayBuffers> ManagedQuery::results() {
    if (is_empty_query()) {
        query_submitted_ = false;
        return buffers_;
    }

    if (!query_submitted_) {
        throw TileDBSOMAError(fmt::format(
            "[ManagedQuery][{}] submit query before reading results", name_));
    }
    query_submitted_ = false;

    // Poll status until query is not INPROGRESS
    Query::Status status;
    do {
        status = query_->query_status();
    } while (status == Query::Status::INPROGRESS);

    LOG_DEBUG(fmt::format(
        "[ManagedQuery] [{}] Query status = {}", name_, (int)status));

    if (status == Query::Status::FAILED) {
        throw TileDBSOMAError(
            fmt::format("[ManagedQuery] [{}] Query FAILED", name_));
    }

    // If the query was ever incomplete, the result buffers contents are not
    // complete.
    if (status == Query::Status::INCOMPLETE) {
        results_complete_ = false;
    }

    // Update ColumnBuffer size to match query results
    size_t num_cells = 0;
    for (auto& name : buffers_->names()) {
        num_cells = buffers_->at(name)->update_size(*query_);
        LOG_DEBUG(fmt::format(
            "[ManagedQuery] [{}] Buffer {} cells={}", name_, name, num_cells));
    }
    total_num_cells_ += num_cells;

    // TODO: retry the query with larger buffers
    if (status == Query::Status::INCOMPLETE && !num_cells) {
        throw TileDBSOMAError(
            fmt::format("[ManagedQuery] [{}] Buffers are too small.", name_));
    }

    return buffers_;
}

};  // namespace tiledbsoma
