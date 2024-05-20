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

#include "managed_query.h"
#include <tiledb/array_experimental.h>
#include <tiledb/attribute_experimental.h>
#include "../utils/logger.h"
#include "utils/common.h"
namespace tiledbsoma {

using namespace tiledb;

//===================================================================
//= public non-static
//===================================================================

ManagedQuery::ManagedQuery(
    std::shared_ptr<Array> array,
    std::shared_ptr<Context> ctx,
    std::string_view name)
    : array_(array)
    , ctx_(ctx)
    , name_(name)
    , schema_(std::make_shared<ArraySchema>(array->schema())) {
    reset();
}

void ManagedQuery::close() {
    if (query_future_.valid()) {
        query_future_.get();
    }
    array_->close();
}

void ManagedQuery::reset() {
    query_ = std::make_unique<Query>(*ctx_, *array_);
    subarray_ = std::make_unique<Subarray>(*ctx_, *array_);

    subarray_range_set_ = false;
    subarray_range_empty_ = {};
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

void ManagedQuery::set_column_data(
    std::shared_ptr<ColumnBuffer> column_buffer) {
    auto column_name = std::string(column_buffer->name());
    bool has_attr = schema_->has_attribute(column_name);
    bool is_sparse = array_->schema().array_type() == TILEDB_SPARSE;

    if (is_sparse) {
        auto data = column_buffer->data<std::byte>();
        query_->set_data_buffer(
            column_name, (void*)data.data(), column_buffer->data_size());
        if (column_buffer->is_var()) {
            auto offsets = column_buffer->offsets();
            query_->set_offsets_buffer(
                column_name, offsets.data(), offsets.size());
        }
        if (column_buffer->is_nullable()) {
            auto validity = column_buffer->validity();
            query_->set_validity_buffer(
                column_name, validity.data(), validity.size());
        }
    } else {
        if (has_attr) {
            auto data = column_buffer->data<std::byte>();
            query_->set_data_buffer(
                column_name, (void*)data.data(), column_buffer->data_size());
            if (column_buffer->is_var()) {
                auto offsets = column_buffer->offsets();
                query_->set_offsets_buffer(
                    column_name, offsets.data(), offsets.size());
            }
            if (column_buffer->is_nullable()) {
                auto validity = column_buffer->validity();
                query_->set_validity_buffer(
                    column_name, validity.data(), validity.size());
            }
        } else {
            switch (column_buffer->type()) {
                case TILEDB_STRING_ASCII:
                case TILEDB_STRING_UTF8:
                case TILEDB_CHAR:
                case TILEDB_BLOB:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<std::string>()[0],
                        column_buffer->data<std::string>()[1]);
                    break;
                case TILEDB_FLOAT32:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<float>()[0],
                        column_buffer->data<float>()[1]);
                    break;
                case TILEDB_FLOAT64:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<double>()[0],
                        column_buffer->data<double>()[1]);
                    break;
                case TILEDB_UINT8:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<uint8_t>()[0],
                        column_buffer->data<uint8_t>()[1]);
                    break;
                case TILEDB_INT8:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<int8_t>()[0],
                        column_buffer->data<int8_t>()[1]);
                    break;
                case TILEDB_UINT16:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<uint16_t>()[0],
                        column_buffer->data<uint16_t>()[1]);
                    break;
                case TILEDB_INT16:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<int16_t>()[0],
                        column_buffer->data<int16_t>()[1]);
                    break;
                case TILEDB_UINT32:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<uint32_t>()[0],
                        column_buffer->data<uint32_t>()[1]);
                    break;
                case TILEDB_INT32:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<int32_t>()[0],
                        column_buffer->data<int32_t>()[1]);
                    break;
                case TILEDB_UINT64:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<uint64_t>()[0],
                        column_buffer->data<uint64_t>()[1]);
                    break;
                case TILEDB_INT64:
                case TILEDB_TIME_SEC:
                case TILEDB_TIME_MS:
                case TILEDB_TIME_US:
                case TILEDB_TIME_NS:
                case TILEDB_DATETIME_SEC:
                case TILEDB_DATETIME_MS:
                case TILEDB_DATETIME_US:
                case TILEDB_DATETIME_NS:
                    subarray_->add_range(
                        column_name,
                        column_buffer->data<int64_t>()[0],
                        column_buffer->data<int64_t>()[1]);
                    break;
                default:
                    break;
            }
            query_->set_subarray(*subarray_);
        }
    }
}

void ManagedQuery::setup_read() {
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
        if (array_->schema().array_type() == TILEDB_SPARSE) {
            for (const auto& dim : array_->schema().domain().dimensions()) {
                columns_.push_back(dim.name());
            }
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
}

void ManagedQuery::submit_write(bool sort_coords) {
    if (array_->schema().array_type() == TILEDB_DENSE) {
        query_->set_subarray(*subarray_);
    } else {
        query_->set_layout(
            sort_coords ? TILEDB_UNORDERED : TILEDB_GLOBAL_ORDER);
    }

    if (query_->query_layout() == TILEDB_GLOBAL_ORDER) {
        query_->submit_and_finalize();
    } else {
        query_->submit();
        query_->finalize();
    }
}

void ManagedQuery::submit_read() {
    query_submitted_ = true;
    query_future_ = std::async(std::launch::async, [&]() {
        LOG_DEBUG("[ManagedQuery] submit thread start");
        query_->submit();
        LOG_DEBUG("[ManagedQuery] submit thread done");
    });
}

std::shared_ptr<ArrayBuffers> ManagedQuery::results() {
    if (is_empty_query()) {
        return buffers_;
    }

    if (query_future_.valid()) {
        LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Waiting for query", name_));
        query_future_.wait();
    } else {
        throw TileDBSOMAError(
            fmt::format("[ManagedQuery] [{}] 'query_future_' invalid", name_));
    }

    auto status = query_->query_status();

    if (status == Query::Status::FAILED) {
        throw TileDBSOMAError(
            fmt::format("[ManagedQuery] [{}] Query FAILED", name_));
    }

    // If the query was ever incomplete, the result buffers contents are not
    // complete.
    if (status == Query::Status::INCOMPLETE) {
        results_complete_ = false;
    } else if (status == Query::Status::COMPLETE) {
        results_complete_ = true;
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

    // Visit all attributes and retrieve enumeration vectors
    auto attribute_map = schema_->attributes();
    for (auto& nmit : attribute_map) {
        auto attrname = nmit.first;
        auto attribute = nmit.second;
        auto enumname = AttributeExperimental::get_enumeration_name(
            *ctx_, attribute);
        if (enumname != std::nullopt) {
            auto enumeration = ArrayExperimental::get_enumeration(
                *ctx_, *array_, enumname.value());
            auto enumvec = enumeration.as_vector<std::string>();
            if (!buffers_->contains(attrname)) {
                continue;
            }
            auto colbuf = buffers_->at(attrname);
            colbuf->add_enumeration(enumvec);
            LOG_DEBUG(fmt::format(
                "[ManagedQuery] got Enumeration '{}' for attribute '{}'",
                enumname.value(),
                attrname));
        }
    }
    return buffers_;
}

void ManagedQuery::check_column_name(const std::string& name) {
    if (!buffers_->contains(name)) {
        throw TileDBSOMAError(fmt::format(
            "[ManagedQuery] Column '{}' is not available in the query "
            "results.",
            name));
    }
}
};  // namespace tiledbsoma
