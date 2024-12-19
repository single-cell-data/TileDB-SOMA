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
#include "soma_array.h"
#include "utils/common.h"
#include "utils/logger.h"
#include "utils/util.h"
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

ManagedQuery::ManagedQuery(
    std::unique_ptr<SOMAArray> array,
    std::shared_ptr<Context> ctx,
    std::string_view name)
    : array_(array->arr_)
    , ctx_(ctx)
    , name_(name)
    , schema_(std::make_shared<ArraySchema>(array->arr_->schema())) {
    reset();
}

void ManagedQuery::close() {
    array_->close();
}

void ManagedQuery::reset() {
    query_ = std::make_unique<Query>(*ctx_, *array_);
    subarray_ = std::make_unique<Subarray>(*ctx_, *array_);

    subarray_range_set_ = {};
    subarray_range_empty_ = {};
    columns_.clear();
    results_complete_ = true;
    total_num_cells_ = 0;
    buffers_.reset();
    query_submitted_ = false;
}

void ManagedQuery::set_layout(ResultOrder layout) {
    switch (layout) {
        case ResultOrder::automatic:
            if (schema_->array_type() == TILEDB_SPARSE)
                query_->set_layout(TILEDB_UNORDERED);
            else
                query_->set_layout(TILEDB_ROW_MAJOR);
            break;
        case ResultOrder::rowmajor:
            query_->set_layout(TILEDB_ROW_MAJOR);
            break;
        case ResultOrder::colmajor:
            query_->set_layout(TILEDB_COL_MAJOR);
            break;
        default:
            throw std::invalid_argument(fmt::format(
                "[ManagedQuery] invalid ResultOrder({}) passed",
                static_cast<int>(layout)));
    }
}

void ManagedQuery::select_columns(
    const std::vector<std::string>& names, bool if_not_empty, bool replace) {
    // Return if we are selecting all columns (columns_ is empty) and we want to
    // continue selecting all columns (if_not_empty == true).
    if (if_not_empty && columns_.empty()) {
        return;
    }

    if (replace) {
        reset_columns();
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

void ManagedQuery::reset_columns() {
    columns_.clear();
}

void ManagedQuery::setup_read() {
    // If the query is complete, return so we do not submit it again
    auto status = query_->query_status();
    if (status == Query::Status::COMPLETE) {
        return;
    }

    auto schema = array_->schema();

    _fill_in_subarrays_if_dense(true);

    // If the query is uninitialized, set the subarray for the query
    if (status == Query::Status::UNINITIALIZED) {
        // Set the subarray for range slicing
        query_->set_subarray(*subarray_);
    }

    // If no columns were selected, select all columns.
    // Add dims and attrs in the same order as specified in the schema
    if (columns_.empty()) {
        if (schema.array_type() == TILEDB_SPARSE) {
            for (const auto& dim : schema.domain().dimensions()) {
                columns_.push_back(dim.name());
            }
        }
        int attribute_num = schema.attribute_num();
        for (int i = 0; i < attribute_num; i++) {
            columns_.push_back(schema.attribute(i).name());
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
    _fill_in_subarrays_if_dense(false);

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
        try {
            query_->submit();
        } catch (const std::exception& e) {
            return StatusAndException(false, e.what());
        }
        LOG_DEBUG("[ManagedQuery] submit thread done");
        return StatusAndException(true, "success");
    });
}

std::optional<std::shared_ptr<ArrayBuffers>> ManagedQuery::read_next() {
    setup_read();

    if (is_empty_query() && !query_submitted_) {
        query_submitted_ = true;
        return buffers_;
    }

    if (is_complete(false)) {
        return std::nullopt;
    }

    query_submitted_ = true;
    query_future_ = std::async(std::launch::async, [&]() {
        LOG_DEBUG("[ManagedQuery] submit thread start");
        try {
            query_->submit();
        } catch (const std::exception& e) {
            return StatusAndException(false, e.what());
        }
        LOG_DEBUG("[ManagedQuery] submit thread done");
        return StatusAndException(true, "success");
    });

    if (query_future_.valid()) {
        LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Waiting for query", name_));
        query_future_.wait();
        LOG_DEBUG(
            fmt::format("[ManagedQuery] [{}] Done waiting for query", name_));

        auto retval = query_future_.get();
        if (!retval.succeeded()) {
            throw TileDBSOMAError(fmt::format(
                "[ManagedQuery] [{}] Query FAILED: {}",
                name_,
                retval.message()));
        }

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

// Please see the header-file comments for context.
void ManagedQuery::_fill_in_subarrays_if_dense(bool is_read) {
    LOG_TRACE("[ManagedQuery] _fill_in_subarrays enter");
    // Don't do this on next-page etc.
    if (query_->query_status() != Query::Status::UNINITIALIZED) {
        LOG_TRACE("[ManagedQuery] _fill_in_subarrays exit: initialized");
        return;
    }
    auto schema = array_->schema();

    // Do this only for dense arrays.
    if (schema.array_type() != TILEDB_DENSE) {
        LOG_TRACE("[ManagedQuery] _fill_in_subarrays exit: non-dense");
        return;
    }

    // Don't do this if the array doesn't have new shape AKA current domain.
    auto current_domain = tiledb::ArraySchemaExperimental::current_domain(
        *ctx_, schema);
    if (current_domain.is_empty()) {
        _fill_in_subarrays_if_dense_without_new_shape(is_read);
    } else {
        _fill_in_subarrays_if_dense_with_new_shape(current_domain, is_read);
    }
    LOG_TRACE("[ManagedQuery] _fill_in_subarrays exit");
}

void ManagedQuery::_fill_in_subarrays_if_dense_without_new_shape(bool is_read) {
    LOG_TRACE(
        "[ManagedQuery] _fill_in_subarrays_if_dense_without_new_shape enter");
    // Dense array must have a subarray set for read. If the array is dense and
    // no ranges have been set, add a range for the array's entire non-empty
    // domain on dimension 0. In the case that the non-empty domain does not
    // exist (when the array has not been written to yet), use dimension 0's
    // full domain
    if (_has_any_subarray_range_set()) {
        return;
    }

    std::pair<int64_t, int64_t> array_shape;
    auto schema = array_->schema();

    if (!is_read && subarray_->range_num(0) > 0) {
        LOG_TRACE(
            "[ManagedQuery] _fill_in_subarrays_if_dense_without_new_shape "
            "range 0 is set");
        return;
    }

    if (is_read) {
        // Check if the array has been written to by using the C API as
        // there is no way to to check for an empty domain using the current
        // C++ API.
        int32_t is_empty;
        int64_t ned[2];
        ctx_->handle_error(tiledb_array_get_non_empty_domain_from_index(
            ctx_->ptr().get(), array_->ptr().get(), 0, &ned, &is_empty));

        if (is_empty == 1) {
            array_shape = schema.domain().dimension(0).domain<int64_t>();
        } else {
            array_shape = std::make_pair(ned[0], ned[1]);
        }
    } else {
        // Non-empty domain is not avaiable for access at write time.
        array_shape = schema.domain().dimension(0).domain<int64_t>();
    }

    subarray_->add_range(0, array_shape.first, array_shape.second);
    LOG_TRACE(fmt::format(
        "[ManagedQuery] Add full range to dense subarray dim0 = ({}, {})",
        array_shape.first,
        array_shape.second));

    // Set the subarray for range slicing
    query_->set_subarray(*subarray_);
}

void ManagedQuery::_fill_in_subarrays_if_dense_with_new_shape(
    const CurrentDomain& current_domain, bool is_read) {
    LOG_TRACE(
        "[ManagedQuery] _fill_in_subarrays_if_dense_with_new_shape enter");
    if (current_domain.type() != TILEDB_NDRECTANGLE) {
        throw TileDBSOMAError("found non-rectangle current-domain type");
    }
    NDRectangle ndrect = current_domain.ndrectangle();

    // Loop over dims and apply subarray ranges if not already done by the
    // caller.
    auto schema = array_->schema();
    for (const auto& dim : schema.domain().dimensions()) {
        std::string dim_name = dim.name();
        if (subarray_range_set_[dim_name]) {
            LOG_TRACE(fmt::format(
                "[ManagedQuery] _fill_in_subarrays continue {}", dim_name));
            continue;
        }

        // Dense arrays are (as of this writing in 1.15.0) all DenseNDArray.
        // Per the spec DenseNDArray must only have dims named
        // soma_dim_{i} with i=0,1,2,...,n-1, of type int64.
        if (dim_name.rfind("soma_dim_", 0) != 0) {
            throw TileDBSOMAError(fmt::format(
                "found dense array with unexpected dim name {}", dim_name));
        }
        if (dim.type() != TILEDB_INT64) {
            throw TileDBSOMAError(fmt::format(
                "expected dense arrays to have int64 dims; got {} for {}",
                tiledb::impl::to_str(dim.type()),
                dim_name));
        }

        std::array<int64_t, 2> lo_hi_arr = ndrect.range<int64_t>(dim_name);
        int64_t cd_lo = lo_hi_arr[0];
        int64_t cd_hi = lo_hi_arr[1];
        int64_t lo = cd_lo;
        int64_t hi = cd_hi;

        LOG_TRACE(fmt::format(
            "[ManagedQuery] _fill_in_subarrays_if_dense_with_new_shape dim "
            "name {} current domain ({}, {})",
            dim_name,
            cd_lo,
            cd_hi));

        if (is_read) {
            // Check if the array has been written to by using the C API as
            // there is no way to to check for an empty domain using the current
            // C++ API.
            //
            // Non-empty domain is not avaliable at write time (the core array
            // isn't open for read) -- but, we don't need to do this calculation
            // for writes anyway.
            int32_t is_empty;
            int64_t ned[2];
            ctx_->handle_error(tiledb_array_get_non_empty_domain_from_name(
                ctx_->ptr().get(),
                array_->ptr().get(),
                dim_name.c_str(),
                &ned,
                &is_empty));
            if (is_empty == 1) {
                LOG_TRACE(fmt::format(
                    "[ManagedQuery] _fill_in_subarrays_if_dense_with_new_shape "
                    "dim name {} non-empty domain is absent",
                    dim_name));
            } else {
                int64_t ned_lo = ned[0];
                int64_t ned_hi = ned[1];
                LOG_TRACE(fmt::format(
                    "[ManagedQuery] _fill_in_subarrays_if_dense_with_new_shape "
                    "dim name {} non-empty domain ({}, {})",
                    dim_name,
                    ned_lo,
                    ned_hi));
                if (ned_lo > cd_lo) {
                    lo = ned_lo;
                }
                if (ned_hi < cd_hi) {
                    hi = ned_hi;
                }
            }
        }

        LOG_TRACE(fmt::format(
            "[ManagedQuery] _fill_in_subarrays_if_dense_with_new_shape dim "
            "name {} select ({}, {})",
            dim_name,
            lo,
            hi));

        std::pair<int64_t, int64_t> lo_hi_pair(lo, hi);
        select_ranges(dim_name, std::vector({lo_hi_pair}));
    }
}

std::shared_ptr<ArrayBuffers> ManagedQuery::results() {
    if (is_empty_query()) {
        query_submitted_ = true;
        return buffers_;
    }

    if (query_future_.valid()) {
        LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Waiting for query", name_));
        query_future_.wait();
        LOG_DEBUG(
            fmt::format("[ManagedQuery] [{}] Done waiting for query", name_));

        auto retval = query_future_.get();
        if (!retval.succeeded()) {
            throw TileDBSOMAError(fmt::format(
                "[ManagedQuery] [{}] Query FAILED: {}",
                name_,
                retval.message()));
        }

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

uint64_t ManagedQuery::_get_max_capacity(tiledb_datatype_t index_type) {
    switch (index_type) {
        case TILEDB_INT8:
            return std::numeric_limits<int8_t>::max();
        case TILEDB_UINT8:
            return std::numeric_limits<uint8_t>::max();
        case TILEDB_INT16:
            return std::numeric_limits<int16_t>::max();
        case TILEDB_UINT16:
            return std::numeric_limits<uint16_t>::max();
        case TILEDB_INT32:
            return std::numeric_limits<int32_t>::max();
        case TILEDB_UINT32:
            return std::numeric_limits<uint32_t>::max();
        case TILEDB_INT64:
            return std::numeric_limits<int64_t>::max();
        case TILEDB_UINT64:
            return std::numeric_limits<uint64_t>::max();
        default:
            throw TileDBSOMAError(
                "Saw invalid enumeration index type when trying to extend "
                "enumeration");
    }
}

ArraySchemaEvolution ManagedQuery::_make_se() {
    ArraySchemaEvolution se(*ctx_);
    return se;
}

void ManagedQuery::set_array_data(
    std::unique_ptr<ArrowSchema> arrow_schema,
    std::unique_ptr<ArrowArray> arrow_array) {
    // Go through all columns in the ArrowTable and cast the values to what is
    // in the ArraySchema on disk
    ArraySchemaEvolution se = _make_se();
    bool evolve_schema = false;
    for (auto i = 0; i < arrow_schema->n_children; ++i) {
        bool enmr_extended = _cast_column(
            arrow_schema->children[i], arrow_array->children[i], se);
        evolve_schema = evolve_schema || enmr_extended;
    }
    if (evolve_schema) {
        se.array_evolve(array_->uri());
    }
};

bool ManagedQuery::_cast_column(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se) {
    auto user_type = ArrowAdapter::to_tiledb_format(schema->format);
    bool has_attr = schema_->has_attribute(schema->name);

    // If the attribute is enumerated, but the provided column is not, error out
    if (has_attr && attr_has_enum(schema->name)) {
        if (schema->dictionary == nullptr || array->dictionary == nullptr) {
            throw std::invalid_argument(
                "[SOMAArray] " + std::string(schema->name) +
                " requires dictionary entry");
        }
    }

    // If the attribute is not enumerated, but the provided column is, then we
    // need to use the dictionary values when writing to the array
    if (has_attr && !attr_has_enum(schema->name)) {
        if (schema->dictionary != nullptr && array->dictionary != nullptr) {
            _promote_indexes_to_values(schema, array);

            // Return false because we do not extend the enumeration
            return false;
        }
    }

    // In the general cases that do not apply to the two cases above, we need to
    // cast the passed-in column to be what is the type in the schema on disk.
    // Here we identify the passed-in column type (UserType).
    //
    // If _cast_column_aux extended the enumeration then return true. Otherwise,
    // false
    switch (user_type) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            return _cast_column_aux<std::string>(schema, array, se);
        case TILEDB_BOOL:
            return _cast_column_aux<bool>(schema, array, se);
        case TILEDB_INT8:
            return _cast_column_aux<int8_t>(schema, array, se);
        case TILEDB_UINT8:
            return _cast_column_aux<uint8_t>(schema, array, se);
        case TILEDB_INT16:
            return _cast_column_aux<int16_t>(schema, array, se);
        case TILEDB_UINT16:
            return _cast_column_aux<uint16_t>(schema, array, se);
        case TILEDB_INT32:
            return _cast_column_aux<int32_t>(schema, array, se);
        case TILEDB_UINT32:
            return _cast_column_aux<uint32_t>(schema, array, se);
        case TILEDB_INT64:
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_TIME_HR:
        case TILEDB_TIME_MIN:
        case TILEDB_TIME_SEC:
        case TILEDB_TIME_MS:
        case TILEDB_TIME_US:
        case TILEDB_TIME_NS:
        case TILEDB_TIME_PS:
        case TILEDB_TIME_FS:
        case TILEDB_TIME_AS:
            return _cast_column_aux<int64_t>(schema, array, se);
        case TILEDB_UINT64:
            return _cast_column_aux<uint64_t>(schema, array, se);
        case TILEDB_FLOAT32:
            return _cast_column_aux<float>(schema, array, se);
        case TILEDB_FLOAT64:
            return _cast_column_aux<double>(schema, array, se);
        default:
            throw TileDBSOMAError(fmt::format(
                "Saw invalid TileDB user type when attempting to cast table: "
                "{}",
                tiledb::impl::type_to_str(user_type)));
    }
}

void ManagedQuery::_promote_indexes_to_values(
    ArrowSchema* schema, ArrowArray* array) {
    // This is a column with a dictionary. However, the associated TileDB
    // attribute on disk is not enumerated. We will need to map the dictionary
    // indexes to the associated dictionary values and write the values to disk.
    // Here, we identify the passed-in column type

    auto value_type = ArrowAdapter::to_tiledb_format(
        schema->dictionary->format);
    switch (value_type) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            return _cast_dictionary_values<std::string>(schema, array);
        case TILEDB_BOOL:
            return _cast_dictionary_values<bool>(schema, array);
        case TILEDB_INT8:
            return _cast_dictionary_values<int8_t>(schema, array);
        case TILEDB_UINT8:
            return _cast_dictionary_values<uint8_t>(schema, array);
        case TILEDB_INT16:
            return _cast_dictionary_values<int16_t>(schema, array);
        case TILEDB_UINT16:
            return _cast_dictionary_values<uint16_t>(schema, array);
        case TILEDB_INT32:
            return _cast_dictionary_values<int32_t>(schema, array);
        case TILEDB_UINT32:
            return _cast_dictionary_values<uint32_t>(schema, array);
        case TILEDB_INT64:
            return _cast_dictionary_values<int64_t>(schema, array);
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_TIME_HR:
        case TILEDB_TIME_MIN:
        case TILEDB_TIME_SEC:
        case TILEDB_TIME_MS:
        case TILEDB_TIME_US:
        case TILEDB_TIME_NS:
        case TILEDB_TIME_PS:
        case TILEDB_TIME_FS:
        case TILEDB_TIME_AS:
        case TILEDB_UINT64:
            return _cast_dictionary_values<uint64_t>(schema, array);
        case TILEDB_FLOAT32:
            return _cast_dictionary_values<float>(schema, array);
        case TILEDB_FLOAT64:
            return _cast_dictionary_values<double>(schema, array);
        default:
            throw TileDBSOMAError(fmt::format(
                "Saw invalid TileDB value type when attempting to promote "
                "indexes to values: {}",
                tiledb::impl::type_to_str(value_type)));
    }
}

template <typename T>
void ManagedQuery::_cast_dictionary_values(
    ArrowSchema* schema, ArrowArray* array) {
    // This is a column with a dictionary. However, the associated TileDB
    // attribute on disk is not enumerated. Here, we map the dictionary indexes
    // to the associated dictionary values and set the buffers to use the
    // dictionary values to write to disk. Note the specialized templates for
    // string and Boolean types below

    auto value_array = array->dictionary;

    T* valbuf;
    if (value_array->n_buffers == 3) {
        valbuf = (T*)value_array->buffers[2];
    } else {
        valbuf = (T*)value_array->buffers[1];
    }
    std::vector<T> values(valbuf, valbuf + value_array->length);

    std::vector<int64_t> indexes = _get_index_vector(schema, array);

    std::vector<T> index_to_value;
    for (auto i : indexes) {
        index_to_value.push_back(values[i]);
    }

    setup_write_column(
        schema->name,
        array->length,
        (const void*)index_to_value.data(),
        (uint64_t*)nullptr,
        (uint8_t*)value_array->buffers[0]);
}

template <>
void ManagedQuery::_cast_dictionary_values<std::string>(
    ArrowSchema* schema, ArrowArray* array) {
    // String types require special handling due to large vs regular
    // string/binary

    auto value_schema = schema->dictionary;
    auto value_array = array->dictionary;

    uint64_t num_elems = value_array->length;

    std::vector<uint64_t> offsets_v;
    if ((strcmp(value_schema->format, "U") == 0) ||
        (strcmp(value_schema->format, "Z") == 0)) {
        uint64_t* offsets = (uint64_t*)value_array->buffers[1];
        offsets_v.resize(num_elems + 1);
        offsets_v.assign(offsets, offsets + num_elems + 1);
    } else {
        uint32_t* offsets = (uint32_t*)value_array->buffers[1];
        std::vector<uint32_t> offset_holder(offsets, offsets + num_elems + 1);
        for (auto offset : offset_holder) {
            offsets_v.push_back((uint64_t)offset);
        }
    }

    char* data = (char*)value_array->buffers[2];
    std::string data_v(data, data + offsets_v[offsets_v.size() - 1]);

    std::vector<std::string> values;
    for (size_t offset_idx = 0; offset_idx < offsets_v.size() - 1;
         ++offset_idx) {
        auto beg = offsets_v[offset_idx];
        auto sz = offsets_v[offset_idx + 1] - beg;
        values.push_back(data_v.substr(beg, sz));
    }

    std::vector<int64_t> indexes = ManagedQuery::_get_index_vector(
        schema, array);

    uint64_t offset_sum = 0;
    std::vector<uint64_t> value_offsets = {0};
    std::string index_to_value;
    for (auto i : indexes) {
        auto value = values[i];
        offset_sum += value.size();
        value_offsets.push_back(offset_sum);
        index_to_value.insert(index_to_value.end(), value.begin(), value.end());
    }

    setup_write_column(
        schema->name,
        value_offsets.size() - 1,
        (const void*)index_to_value.data(),
        (uint64_t*)value_offsets.data(),
        (uint8_t*)value_array->buffers[0]);
}

template <>
void ManagedQuery::_cast_dictionary_values<bool>(
    ArrowSchema* schema, ArrowArray* array) {
    // Boolean types require special handling due to bit vs uint8_t
    // representation in Arrow vs TileDB respectively

    auto value_schema = schema->dictionary;
    auto value_array = array->dictionary;

    std::vector<int64_t> indexes = _get_index_vector(schema, array);
    std::vector<uint8_t> values = util::cast_bit_to_uint8(
        value_schema, value_array);
    std::vector<uint8_t> index_to_value;

    for (auto i : indexes) {
        index_to_value.push_back(values[i]);
    }

    setup_write_column(
        schema->name,
        array->length,
        (const void*)index_to_value.data(),
        (uint64_t*)nullptr,
        (uint8_t*)value_array->buffers[0]);
}

template <typename UserType>
bool ManagedQuery::_cast_column_aux(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se) {
    // We need to cast the passed-in column to be what is the type in the schema
    // on disk. Here we identify the on-disk attribute or dimension type
    // (DiskType).

    tiledb_datatype_t disk_type;
    std::string name(schema->name);
    if (schema_->has_attribute(name)) {
        disk_type = schema_->attribute(name).type();
    } else {
        disk_type = schema_->domain().dimension(name).type();
    }

    // If _set_column extended the enumeration then return true. Otherwise,
    // false
    switch (disk_type) {
        case TILEDB_BOOL:
        case TILEDB_INT8:
            return _set_column<UserType, int8_t>(schema, array, se);
        case TILEDB_UINT8:
            return _set_column<UserType, uint8_t>(schema, array, se);
        case TILEDB_INT16:
            return _set_column<UserType, int16_t>(schema, array, se);
        case TILEDB_UINT16:
            return _set_column<UserType, uint16_t>(schema, array, se);
        case TILEDB_INT32:
            return _set_column<UserType, int32_t>(schema, array, se);
        case TILEDB_UINT32:
            return _set_column<UserType, uint32_t>(schema, array, se);
        case TILEDB_INT64:
        case TILEDB_DATETIME_YEAR:
        case TILEDB_DATETIME_MONTH:
        case TILEDB_DATETIME_WEEK:
        case TILEDB_DATETIME_DAY:
        case TILEDB_DATETIME_HR:
        case TILEDB_DATETIME_MIN:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_DATETIME_PS:
        case TILEDB_DATETIME_FS:
        case TILEDB_DATETIME_AS:
        case TILEDB_TIME_HR:
        case TILEDB_TIME_MIN:
        case TILEDB_TIME_SEC:
        case TILEDB_TIME_MS:
        case TILEDB_TIME_US:
        case TILEDB_TIME_NS:
        case TILEDB_TIME_PS:
        case TILEDB_TIME_FS:
        case TILEDB_TIME_AS:
            return _set_column<UserType, int64_t>(schema, array, se);
        case TILEDB_UINT64:
            return _set_column<UserType, uint64_t>(schema, array, se);
        case TILEDB_FLOAT32:
            return _set_column<UserType, float>(schema, array, se);
        case TILEDB_FLOAT64:
            return _set_column<UserType, double>(schema, array, se);
        default:
            throw TileDBSOMAError(
                "Saw invalid TileDB disk type when attempting to cast "
                "column: " +
                tiledb::impl::type_to_str(disk_type));
    }
}

template <>
bool ManagedQuery::_cast_column_aux<std::string>(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se) {
    (void)se;  // se is unused in std::string specialization

    // A few things in play here:
    // * Whether the column (array) has 3 buffers (validity, offset, data)
    //   or 2 (validity, data).
    // * The data is always char* and the validity is always uint8*
    //   but the offsets are 32-bit or 64-bit.
    // * The array-level offset might not be zero. (This happens
    //   when people pass of things like arrow_table[:n] or arrow_table[n:]
    //   from Python/R.)

    if (array->n_buffers != 3) {
        throw TileDBSOMAError(fmt::format(
            "[ManagedQuery] internal error: Arrow-table string column should "
            "have 3 buffers; got {}",
            array->n_buffers));
    }

    const char* data = (const char*)array->buffers[2];
    uint8_t* validity = (uint8_t*)array->buffers[0];

    // If this is a table-slice, slice into the validity buffer.
    if (validity != nullptr) {
        validity += array->offset;
    }
    // If this is a table-slice, do *not* slice into the data
    // buffer since it is indexed via offsets, which we slice
    // into below.

    if ((strcmp(schema->format, "U") == 0) ||
        (strcmp(schema->format, "Z") == 0)) {
        // If this is a table-slice, slice into the offsets buffer.
        uint64_t* offset = (uint64_t*)array->buffers[1] + array->offset;
        setup_write_column(
            schema->name, array->length, (const void*)data, offset, validity);

    } else {
        // If this is a table-slice, slice into the offsets buffer.
        uint32_t* offset = (uint32_t*)array->buffers[1] + array->offset;
        setup_write_column(
            schema->name, array->length, (const void*)data, offset, validity);
    }
    return false;
}

template <>
bool ManagedQuery::_cast_column_aux<bool>(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se) {
    (void)se;  // se is unused in bool specialization

    auto casted = util::cast_bit_to_uint8(schema, array);
    setup_write_column(
        schema->name,
        array->length,
        (const void*)casted.data(),
        (uint64_t*)nullptr,
        (uint8_t*)array->buffers[0]);
    return false;
}

bool ManagedQuery::_extend_enumeration(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    ArraySchemaEvolution se) {
    // For columns with dictionaries, we need to identify whether the

    auto enmr = ArrayExperimental::get_enumeration(
        *ctx_, *array_, index_schema->name);
    auto value_type = enmr.type();

    switch (value_type) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            return _extend_and_evolve_schema<std::string>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_INT8:
            return _extend_and_evolve_schema<int8_t>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_BOOL:
        case TILEDB_UINT8:
            return _extend_and_evolve_schema<uint8_t>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_INT16:
            return _extend_and_evolve_schema<int16_t>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_UINT16:
            return _extend_and_evolve_schema<uint16_t>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_INT32:
            return _extend_and_evolve_schema<int32_t>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_UINT32:
            return _extend_and_evolve_schema<uint32_t>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_INT64:
            return _extend_and_evolve_schema<int64_t>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_UINT64:
            return _extend_and_evolve_schema<uint64_t>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_FLOAT32:
            return _extend_and_evolve_schema<float>(
                value_schema, value_array, index_schema, index_array, se);
        case TILEDB_FLOAT64:
            return _extend_and_evolve_schema<double>(
                value_schema, value_array, index_schema, index_array, se);
        default:
            throw TileDBSOMAError(fmt::format(
                "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                tiledb::impl::type_to_str(value_type)));
    }
}

template <typename ValueType>
bool ManagedQuery::_extend_and_evolve_schema(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    ArraySchemaEvolution se) {
    // We need to check if we are writing any new enumeration values. If so,
    // extend and evolve the schema. If not, just set the write buffers to the
    // dictionary's indexes as-is

    // Get all the enumeration values in the passed-in column
    std::vector<ValueType> enums_in_write;
    uint64_t num_elems = value_array->length;
    if (strcmp(value_schema->format, "b") == 0) {
        // Specially handle Boolean types as their representation in Arrow (bit)
        // is different from what is in TileDB (uint8_t)
        auto casted = util::cast_bit_to_uint8(value_schema, value_array);
        enums_in_write.assign(
            (ValueType*)casted.data(), (ValueType*)casted.data() + num_elems);
    } else {
        // General case
        const void* data;
        if (value_array->n_buffers == 3) {
            data = value_array->buffers[2];
        } else {
            data = value_array->buffers[1];
        }
        enums_in_write.assign((ValueType*)data, (ValueType*)data + num_elems);
    }

    // Get all the enumeration values in the on-disk TileDB attribute
    std::string column_name = index_schema->name;
    auto enmr = ArrayExperimental::get_enumeration(*ctx_, *array_, column_name);
    std::vector<ValueType> enums_existing = enmr.as_vector<ValueType>();

    // Find any new enumeration values
    std::vector<ValueType> extend_values;
    for (auto enum_val : enums_in_write) {
        if (std::find(enums_existing.begin(), enums_existing.end(), enum_val) ==
            enums_existing.end()) {
            extend_values.push_back(enum_val);
        }
    }

    // extend_values = {true, false};
    if (extend_values.size() != 0) {
        // We have new enumeration values; additional processing needed

        // Check if the number of new enumeration will cause an overflow if
        // extended
        auto disk_index_type = schema_->attribute(column_name).type();
        uint64_t max_capacity = _get_max_capacity(disk_index_type);
        auto free_capacity = max_capacity - enums_existing.size();
        if (free_capacity < extend_values.size()) {
            throw TileDBSOMAError(
                "Cannot extend enumeration; reached maximum capacity");
        }

        // Take the existing enumeration values on disk and extend with the new
        // enumeration values
        auto extended_enmr = enmr.extend(extend_values);
        se.extend_enumeration(extended_enmr);

        // If the passed-in enumerations are only a subset of the new extended
        // enumerations, then we will need to remap the indexes. ie. the user
        // passes in values [B, C] which maps to indexes [0, 1]. However, the
        // full set of extended enumerations is [A, B, C] which means we need to
        // remap [B, C] to be indexes [1, 2]
        ManagedQuery::_remap_indexes(
            column_name,
            extended_enmr,
            enums_in_write,
            index_schema,
            index_array);

        // The enumeration was extended
        return true;
    } else {
        // Example:
        //
        // * Already on storage/schema there are values a,b,c with indices
        //   0,1,2.
        // * User appends values b,c which, within the Arrow data coming in
        //   from the user, have indices 0,1.
        // * We need to remap those to 1,2.
        ManagedQuery::_remap_indexes(
            column_name, enmr, enums_in_write, index_schema, index_array);

        // The enumeration was not extended
        return false;
    }
}

template <>
bool ManagedQuery::_extend_and_evolve_schema<std::string>(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    ArraySchemaEvolution se) {
    uint64_t num_elems = value_array->length;

    std::vector<uint64_t> offsets_v;
    if ((strcmp(value_schema->format, "U") == 0) ||
        (strcmp(value_schema->format, "Z") == 0)) {
        uint64_t* offsets = (uint64_t*)value_array->buffers[1];
        offsets_v.assign(offsets, offsets + num_elems + 1);
    } else {
        uint32_t* offsets = (uint32_t*)value_array->buffers[1];
        for (size_t i = 0; i < num_elems + 1; ++i) {
            offsets_v.push_back((uint64_t)offsets[i]);
        }
    }

    char* data = (char*)value_array->buffers[2];
    std::string data_v(data, data + offsets_v[num_elems]);

    std::vector<std::string> enums_in_write;
    for (size_t i = 0; i < num_elems; ++i) {
        auto beg = offsets_v[i];
        auto sz = offsets_v[i + 1] - beg;
        enums_in_write.push_back(data_v.substr(beg, sz));
    }

    std::string column_name = index_schema->name;
    auto enmr = ArrayExperimental::get_enumeration(*ctx_, *array_, column_name);
    std::vector<std::string> extend_values;
    auto enums_existing = enmr.as_vector<std::string>();
    for (auto enum_val : enums_in_write) {
        if (std::find(enums_existing.begin(), enums_existing.end(), enum_val) ==
            enums_existing.end()) {
            extend_values.push_back(enum_val);
        }
    }

    if (extend_values.size() != 0) {
        // Check that we extend the enumeration values without
        // overflowing
        auto disk_index_type = schema_->attribute(column_name).type();
        uint64_t max_capacity = ManagedQuery::_get_max_capacity(
            disk_index_type);
        auto free_capacity = max_capacity - enums_existing.size();
        if (free_capacity < extend_values.size()) {
            throw TileDBSOMAError(
                "Cannot extend enumeration; reached maximum capacity");
        }

        auto extended_enmr = enmr.extend(extend_values);
        se.extend_enumeration(extended_enmr);

        ManagedQuery::_remap_indexes(
            column_name,
            extended_enmr,
            enums_in_write,
            index_schema,
            index_array);

        return true;
    } else {
        // Example:
        //
        // * Already on storage/schema there are values a,b,c with indices
        //   0,1,2.
        // * User appends values b,c which, within the Arrow data coming in
        //   from the user, have indices 0,1.
        // * We need to remap those to 1,2.

        ManagedQuery::_remap_indexes(
            column_name, enmr, enums_in_write, index_schema, index_array);
    }
    return false;
}
};  // namespace tiledbsoma