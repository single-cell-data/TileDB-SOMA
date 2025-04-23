/**
 * @file   managed_query.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the performing TileDB queries.
 */

#include "managed_query.h"

#include <tiledb/array_experimental.h>
#include <tiledb/attribute_experimental.h>

#include <unordered_set>

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
    : ctx_(ctx)
    , array_(array)
    , name_(name)
    , schema_(std::make_shared<ArraySchema>(array->schema())) {
    reset();
}

ManagedQuery::ManagedQuery(
    SOMAArray array, std::shared_ptr<Context> ctx, std::string_view name)
    : ctx_(ctx)
    , array_(array.arr_)
    , name_(name)
    , schema_(std::make_shared<ArraySchema>(array.arr_->schema())) {
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
        case ResultOrder::unordered:
            query_->set_layout(TILEDB_UNORDERED);
            break;
        case ResultOrder::global:
            query_->set_layout(TILEDB_GLOBAL_ORDER);
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

    layout_ = layout;
}

ResultOrder ManagedQuery::result_order() {
    return layout_;
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

    if (array_->query_type() != TILEDB_WRITE) {
        throw TileDBSOMAError(
            "[ManagedQuery] write requires array to be opened in write mode");
    }

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

    // When we evolve the schema, the ArraySchema needs to be updated to the
    // latest version so re-open the Array
    array_->close();
    array_->open(TILEDB_WRITE);
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
    // no ranges have been set, add a range for the array's entire full
    // domain on dimension 0
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

    array_shape = schema.domain().dimension(0).domain<int64_t>();
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
            auto [dom_lo, dom_hi] = schema.domain()
                                        .dimension(0)
                                        .domain<int64_t>();
            LOG_TRACE(fmt::format(
                "[ManagedQuery] _fill_in_subarrays_if_dense_with_new_shape "
                "dim name {} non-empty domain ({}, {})",
                dim_name,
                dom_lo,
                dom_hi));
            if (dom_lo > cd_lo) {
                lo = dom_lo;
            }
            if (dom_hi < cd_hi) {
                hi = dom_hi;
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
    ArrowSchema* arrow_schema, ArrowArray* arrow_array) {
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
        std::nullopt);  // validities are set by index column
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

    std::string_view data(
        static_cast<const char*>(value_array->buffers[2]), offsets_v.back());

    std::vector<std::string_view> values;
    for (size_t offset_idx = 0; offset_idx < offsets_v.size() - 1;
         ++offset_idx) {
        auto beg = offsets_v[offset_idx];
        auto sz = offsets_v[offset_idx + 1] - beg;
        values.push_back(data.substr(beg, sz));
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
        std::nullopt);  // validities are set by index column
}

template <>
void ManagedQuery::_cast_dictionary_values<bool>(
    ArrowSchema* schema, ArrowArray* array) {
    // Boolean types require special handling due to bit vs uint8_t
    // representation in Arrow vs TileDB respectively

    auto value_schema = schema->dictionary;
    auto value_array = array->dictionary;

    std::vector<int64_t> indexes = _get_index_vector(schema, array);
    std::vector<uint8_t> values = _bool_data_bits_to_bytes(
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
        std::nullopt);  // validities are set by index column
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

    const void* data = array->buffers[2];
    std::optional<std::vector<uint8_t>> validity = _cast_validity_buffer(array);

    // If this is a table-slice, do *not* slice into the data
    // buffer since it is indexed via offsets, which we slice
    // into below.

    if ((strcmp(schema->format, "U") == 0) ||
        (strcmp(schema->format, "Z") == 0)) {
        // If this is a table-slice, slice into the offsets buffer.
        uint64_t* offset = (uint64_t*)array->buffers[1] + array->offset;
        setup_write_column(schema->name, array->length, data, offset, validity);

    } else {
        // If this is a table-slice, slice into the offsets buffer.
        uint32_t* offset = (uint32_t*)array->buffers[1] + array->offset;
        setup_write_column(schema->name, array->length, data, offset, validity);
    }
    return false;
}

template <>
bool ManagedQuery::_cast_column_aux<bool>(
    ArrowSchema* schema, ArrowArray* array, ArraySchemaEvolution se) {
    (void)se;  // se is unused in bool specialization

    auto casted = _bool_data_bits_to_bytes(schema, array);

    setup_write_column(
        schema->name,
        array->length,
        (const void*)casted.data(),
        (uint64_t*)nullptr,
        _cast_validity_buffer(array));
    return false;
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
            throw TileDBSOMAError(
                "internal coding error: template-specialization failure for "
                "boolean in _cast_column_aux");
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            throw TileDBSOMAError(
                "internal coding error: template-specialization failure for "
                "string in _cast_column_aux");

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

bool ManagedQuery::_extend_and_write_enumeration(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    Enumeration enmr,
    ArraySchemaEvolution& se) {
    // For columns with dictionaries, we need to identify the data type of the
    // enumeration to extend any new enumeration values

    auto value_type = enmr.type();

    switch (value_type) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            return _extend_and_evolve_schema_and_write<std::string>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_INT8:
            return _extend_and_evolve_schema_and_write<int8_t>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_BOOL:
        case TILEDB_UINT8:
            return _extend_and_evolve_schema_and_write<uint8_t>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_INT16:
            return _extend_and_evolve_schema_and_write<int16_t>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_UINT16:
            return _extend_and_evolve_schema_and_write<uint16_t>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_INT32:
            return _extend_and_evolve_schema_and_write<int32_t>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_UINT32:
            return _extend_and_evolve_schema_and_write<uint32_t>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_INT64:
            return _extend_and_evolve_schema_and_write<int64_t>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_UINT64:
            return _extend_and_evolve_schema_and_write<uint64_t>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_FLOAT32:
            return _extend_and_evolve_schema_and_write<float>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        case TILEDB_FLOAT64:
            return _extend_and_evolve_schema_and_write<double>(
                value_schema, value_array, index_schema, index_array, enmr, se);
        default:
            throw TileDBSOMAError(fmt::format(
                "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                tiledb::impl::type_to_str(value_type)));
    }
}

template <>
bool ManagedQuery::_extend_and_evolve_schema_and_write<std::string>(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    Enumeration enmr,
    ArraySchemaEvolution& se) {
    std::string column_name = index_schema->name;

    const auto [was_extended, enum_values_in_write, extended_enmr] =
        _extend_and_evolve_schema_with_details<std::string, std::string_view>(
            value_schema,
            value_array,
            index_schema,
            index_array,
            column_name,
            true,
            enmr,
            se);

    if (was_extended) {
        ManagedQuery::_remap_indexes(
            column_name,
            extended_enmr,
            enum_values_in_write,
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
            column_name, enmr, enum_values_in_write, index_schema, index_array);
    }
    return false;
}

template <typename ValueType>
bool ManagedQuery::_extend_and_evolve_schema_and_write(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    Enumeration enmr,
    ArraySchemaEvolution& se) {
    std::string column_name = index_schema->name;

    const auto [was_extended, enum_values_in_write, extended_enmr] =
        _extend_and_evolve_schema_with_details<ValueType, ValueType>(
            value_schema,
            value_array,
            index_schema,
            index_array,
            column_name,
            true,
            enmr,
            se);

    if (was_extended) {
        // If the passed-in enumerations are only a subset of the new extended
        // enumerations, then we will need to remap the indexes. ie. the user
        // passes in values [B, C] which maps to indexes [0, 1]. However, the
        // full set of extended enumerations is [A, B, C] which means we need to
        // remap [B, C] to be indexes [1, 2]
        ManagedQuery::_remap_indexes(
            column_name,
            extended_enmr,
            enum_values_in_write,
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
            column_name, enmr, enum_values_in_write, index_schema, index_array);

        // The enumeration was not extended
        return false;
    }
}

bool ManagedQuery::_extend_enumeration(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    const std::string& column_name,
    bool deduplicate,
    Enumeration enmr,
    ArraySchemaEvolution& se) {
    // For columns with dictionaries, we need to identify the data type of the
    // enumeration to extend any new enumeration values

    auto value_type_in_schema = enmr.type();
    auto value_type_in_data = ArrowAdapter::to_tiledb_format(
        value_schema->format);

    if (value_type_in_schema != value_type_in_data) {
        throw TileDBSOMAError(fmt::format(
            "extend_enumeration: data type '{}' != schema type '{}'",
            tiledb::impl::type_to_str(value_type_in_data),
            tiledb::impl::type_to_str(value_type_in_schema)));
    }

    switch (value_type_in_schema) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR:
            return _extend_and_evolve_schema_without_details<
                std::string,
                std::string_view>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_INT8:
            return _extend_and_evolve_schema_without_details<int8_t, int8_t>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_BOOL:
        case TILEDB_UINT8:
            return _extend_and_evolve_schema_without_details<uint8_t, uint8_t>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_INT16:
            return _extend_and_evolve_schema_without_details<int16_t, int16_t>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_UINT16:
            return _extend_and_evolve_schema_without_details<
                uint16_t,
                uint16_t>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_INT32:
            return _extend_and_evolve_schema_without_details<int32_t, int32_t>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_UINT32:
            return _extend_and_evolve_schema_without_details<
                uint32_t,
                uint32_t>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_INT64:
            return _extend_and_evolve_schema_without_details<int64_t, int64_t>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_UINT64:
            return _extend_and_evolve_schema_without_details<
                uint64_t,
                uint64_t>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_FLOAT32:
            return _extend_and_evolve_schema_without_details<float, float>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        case TILEDB_FLOAT64:
            return _extend_and_evolve_schema_without_details<double, double>(
                value_schema, value_array, column_name, deduplicate, enmr, se);
        default:
            throw TileDBSOMAError(fmt::format(
                "ArrowAdapter: Unsupported TileDB enumeration datatype: {} ",
                tiledb::impl::type_to_str(value_type_in_schema)));
    }
}

template <>
bool ManagedQuery::
    _extend_and_evolve_schema_without_details<std::string, std::string_view>(
        ArrowSchema* value_schema,
        ArrowArray* value_array,
        const std::string& column_name,
        bool deduplicate,
        Enumeration enmr,
        ArraySchemaEvolution& se) {
    const auto [was_extended, _1, _2] =
        _extend_and_evolve_schema_with_details<std::string, std::string_view>(
            value_schema,
            value_array,
            nullptr,
            nullptr,
            column_name,
            deduplicate,
            enmr,
            se);
    return was_extended;
}

template <typename ValueType, typename ValueViewType>
bool ManagedQuery::_extend_and_evolve_schema_without_details(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    const std::string& column_name,
    bool deduplicate,
    Enumeration enmr,
    ArraySchemaEvolution& se) {
    const auto [was_extended, _1, _2] =
        _extend_and_evolve_schema_with_details<ValueType, ValueType>(
            value_schema,
            value_array,
            nullptr,
            nullptr,
            column_name,
            deduplicate,
            enmr,
            se);
    return was_extended;
}

// We need to check if we are writing any new enumeration values. If so,
// extend and evolve the schema. If not, just set the write buffers to the
// dictionary's indexes as-is
template <>
std::tuple<
    bool,                           // was_extended
    std::vector<std::string_view>,  // enum_values_in_write
    Enumeration>                    // extended_enmr
ManagedQuery::_extend_and_evolve_schema_with_details<std::string>(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,  // null for extend-enum, non-null for write
    ArrowArray* index_array,    // null for extend-enum, non-null for write
    const std::string& column_name,
    bool deduplicate,
    Enumeration enmr,
    ArraySchemaEvolution& se) {
    uint64_t num_elems = value_array->length;

    if (value_array->n_buffers != 3) {
        throw std::invalid_argument(fmt::format(
            "[ManagedQuery] _extend_and_evolve_schema_with_details string: "
            "expected values n_buffers == 3; got {}",
            value_array->n_buffers));
    }

    if (value_array->null_count != 0) {
        throw std::invalid_argument(fmt::format(
            "[ManagedQuery] _extend_and_evolve_schema_with_details string: "
            "null values are not supported"));
    }

    // Set up input values as a char buffer, and offsets within it.  This is
    // zero-copy on the data buffers, since the Arrow data buffers are already
    // contiguous.
    //
    // Arrow var-sized cells can have 32-bit or 64-bit offsets.  TileDB only has
    // 64-bit offsets. Convert from the former to the latter.
    //
    // The + 1 is for the following reason: Suppose the inputs are ["hello",
    // "goodbye"]. Then the first offset is 0, the second is 5, and the third
    // is 12. This makes it possible to locate the end of the string "goodbye"
    // within the char buffer "hellogoodbye". More generally, it takes n+1
    // offsets to specify the starts and ends of n string values within a
    // column.
    std::vector<uint64_t> offsets_v;
    if ((strcmp(value_schema->format, "U") == 0) ||
        (strcmp(value_schema->format, "Z") == 0)) {
        uint64_t* offsets = (uint64_t*)value_array->buffers[1];
        offsets_v.assign(
            offsets + value_array->offset,
            offsets + value_array->offset + num_elems + 1);
    } else {
        uint32_t* offsets = (uint32_t*)value_array->buffers[1];
        for (size_t i = 0; i < num_elems + 1; ++i) {
            offsets_v.push_back((uint64_t)offsets[i + value_array->offset]);
        }
    }

    std::string_view data_as_char(
        static_cast<const char*>(value_array->buffers[2]), offsets_v.back());

    // Create a vector of string-views into the char buffer.
    // We need this in order to partition the requested values
    // into the ones already in the schema, and the ones needing
    // to be added to the schema.
    std::vector<std::string_view> enum_values_in_write;
    std::unordered_set<std::string_view> unique_values_in_write;
    for (size_t i = 0; i < num_elems; ++i) {
        auto beg = offsets_v[i];
        auto sz = offsets_v[i + 1] - beg;
        auto enum_val = data_as_char.substr(beg, sz);
        enum_values_in_write.push_back(enum_val);
        unique_values_in_write.insert(enum_val);
    }

    // Check for non-unique values in the input, even before we check the
    // values in the array schema. See also sc-65078.
    if (enum_values_in_write.size() != unique_values_in_write.size()) {
        // std::range_error maps to Python ValueError
        throw std::range_error(fmt::format(
            "[extend_enumeration] new values provided for column '{}' must "
            "be unique within themselves, irrespective of the deduplicate flag",
            column_name));
    }

    // Separate out the values already in the array schema from the
    // values not already in the array schema.
    auto enum_values_existing = _enumeration_values_view<std::string_view>(
        enmr);
    std::unordered_set<std::string_view> existing_enums_set;
    for (const auto& existing_enum_val : enum_values_existing) {
        existing_enums_set.insert(existing_enum_val);
    }

    std::optional<std::unordered_set<std::string_view>>
        opt_covered_values = _find_covered_enum_values(
            enum_values_in_write, index_schema, index_array);

    std::vector<std::string_view> enum_values_to_add;
    size_t total_size = 0;

    if (opt_covered_values.has_value()) {
        const auto& covered_values = opt_covered_values.value();
        for (size_t i = 0; i < enum_values_in_write.size(); i++) {
            const auto& enum_val = enum_values_in_write[i];
            if (!existing_enums_set.contains(enum_val)) {
                if (covered_values.find(enum_val) != covered_values.end()) {
                    enum_values_to_add.push_back(enum_val);
                    total_size += enum_val.size();
                }
            }
        }
    } else {
        for (size_t i = 0; i < enum_values_in_write.size(); i++) {
            const auto& enum_val = enum_values_in_write[i];
            if (!existing_enums_set.contains(enum_val)) {
                enum_values_to_add.push_back(enum_val);
                total_size += enum_val.size();
            }
        }
    }

    // There are two paths to enumeration extension:
    // 1. The user does dataframe.write.
    //    Here, if the on-schema values are a,b,c and the values being
    //    written are c,d,e, then, the expected UX is that values d,e
    //    will be handed to the core enumeration-extend.
    // 2. The user does dataframe.extend_enumeration_values.
    //    Here, if the on-schema values are a,b,c and the values being
    //    extended are c,d,e, then, the expected UX (requested in sc-63930) is:
    //    * By default, that's an error -- they were supposed to just pass d,e.
    //    * If the deduplicate flag was passed, then we allow that, but
    //      we still need to pass core only the d,e values. (It will
    //      throw otherwise).
    if (!deduplicate &&
        enum_values_to_add.size() != enum_values_in_write.size()) {
        throw TileDBSOMAError(fmt::format(
            "[extend_enumeration] one or more values provided are already "
            "present in the enumeration for column '{}', and deduplicate was "
            "not specified",
            column_name));
    }

    // Extend the enumeration in the schema, if there are any values to be
    // added.
    if (enum_values_to_add.size() != 0) {
        // Check that we extend the enumeration values without
        // overflowing
        auto disk_index_type = schema_->attribute(column_name).type();
        uint64_t max_capacity = ManagedQuery::_get_max_capacity(
            disk_index_type);
        auto free_capacity = max_capacity - enum_values_existing.size();
        if (free_capacity < enum_values_to_add.size()) {
            throw TileDBSOMAError(
                "Cannot extend enumeration; reached maximum capacity");
        }

        std::vector<uint8_t> extend_data(total_size);
        std::vector<uint64_t> extend_offsets;
        extend_offsets.reserve(enum_values_to_add.size());
        size_t curr_offset = 0;
        for (const auto& enum_val : enum_values_to_add) {
            std::memcpy(
                extend_data.data() + curr_offset,
                enum_val.data(),
                enum_val.size());
            extend_offsets.push_back(curr_offset);
            curr_offset += enum_val.size();
        }

        auto extended_enmr = enmr.extend(
            extend_data.data(),
            extend_data.size(),
            extend_offsets.data(),
            enum_values_to_add.size() * sizeof(uint64_t));
        se.extend_enumeration(extended_enmr);

        return std::tuple{
            true,  // was_extended
            enum_values_in_write,
            extended_enmr};

    } else {
        return std::tuple{
            false,  // was_extended
            enum_values_in_write,
            enmr};
    }
}

// We need to check if we are writing any new enumeration values. If so,
// extend and evolve the schema. If not, just set the write buffers to the
// dictionary's indexes as-is
template <typename ValueType, typename ValueViewType>
std::tuple<
    bool,                        // was_extended
    std::vector<ValueViewType>,  // enum_values_in_write
    Enumeration>                 //  extended_enmr
ManagedQuery::_extend_and_evolve_schema_with_details(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,  // null for extend-enum, non-null for write
    ArrowArray* index_array,    // null for extend-enum, non-null for write
    const std::string& column_name,
    bool deduplicate,
    Enumeration enmr,
    ArraySchemaEvolution& se) {
    if (value_array->n_buffers != 2) {
        // Higher-level code should be catching this with an error message which
        // is intended to be user-facing. Here is a low-level check before we
        // dereference buffers[1] and buffers[2].
        throw std::invalid_argument(fmt::format(
            "[ManagedQuery] _extend_and_evolve_schema_with_details non-string: "
            "internal coding error: expected n_buffers == 2; got {}",
            value_array->n_buffers));
    }

    if (value_array->null_count != 0) {
        throw std::invalid_argument(fmt::format(
            "[ManagedQuery] _extend_and_evolve_schema_with_details non-string: "
            "null values are not supported"));
    }
    std::optional<std::vector<uint8_t>> validities = std::nullopt;
    if (index_array != nullptr) {
        validities = _cast_validity_buffer(index_array);
    }

    // Get all the enumeration values in the passed-in column
    std::vector<ValueType> enum_values_in_write;
    uint64_t num_elems = value_array->length;
    if (strcmp(value_schema->format, "b") == 0) {
        // Specially handle Boolean types as their representation in Arrow (bit)
        // is different from what is in TileDB (uint8_t). Here we must copy.
        auto expanded = _bool_data_bits_to_bytes(value_schema, value_array);
        enum_values_in_write.assign(
            expanded.data(), expanded.data() + num_elems);
    } else {
        // General case. This is zero-copy since the Arrow buffer is already
        // contiguous, and the Arrow in-memory storage model is identical to the
        // TileDB in-memory storage model.
        ValueType* data;
        if (value_array->n_buffers == 3) {
            data = (ValueType*)value_array->buffers[2] + value_array->offset;
        } else {
            data = (ValueType*)value_array->buffers[1] + value_array->offset;
        }
        enum_values_in_write.assign(
            (ValueType*)data, (ValueType*)data + num_elems);
    }
    std::unordered_set<ValueType> unique_values_in_write(
        enum_values_in_write.begin(), enum_values_in_write.end());

    // Check for non-unique values in the input, even before we check the
    // values in the array schema. See also sc-65078.
    if (enum_values_in_write.size() != unique_values_in_write.size()) {
        // std::range_error maps to Python ValueError
        throw std::range_error(fmt::format(
            "[extend_enumeration] new values provided for column '{}' must "
            "be unique within themselves, irrespective of the deduplicate flag",
            column_name));
    }

    // The enumeration values are tracked here in as a vector of string-views,
    // even for non-string-valued data.  We do this so that we can
    // differentiate, at the bit-by-bit level, between various floating-point
    // NaN and Inf values.
    size_t n = enum_values_in_write.size();
    std::vector<std::string_view> enum_values_in_write_as_sv(n);
    for (size_t i = 0; i < n; i++) {
        auto sv = std::string_view(
            reinterpret_cast<const char*>(&enum_values_in_write[i]),
            sizeof(enum_values_in_write[i]));
        enum_values_in_write_as_sv[i] = sv;
    }

    // Get all the enumeration values in the on-disk TileDB attribute.
    std::vector<ValueType> enum_values_existing = enmr.as_vector<ValueType>();

    // Separate out the values already in the array schema from the
    // values not already in the array schema.
    //
    // One might think it would be simpler to use
    //
    //   std::unordered_set<ValueType> existing_enums_set;
    //
    // and one would be correct. However, core uses bitwise comparisons, and
    // floating-point NaNs have the following properties: (1) NaN != NaN, and
    // (2) there are multiple floating-point bit patterns which are NaN.  It's
    // simplest to just use the same logic core does, making std::string_view on
    // our elements.
    //
    // Specifically, please see
    // https://github.com/TileDB-Inc/TileDB/blob/2.27.2/tiledb/sm/array_schema/enumeration.cc#L417-L456
    //
    // It is important that we use core's logic here, so that when we are
    // able to access its hashmap directly without constructing our own,
    // that transition will be seamless.
    std::unordered_set<std::string_view> existing_enums_set;
    for (const auto& enum_value_existing : enum_values_existing) {
        auto sv = std::string_view(
            static_cast<char*>((char*)&enum_value_existing),
            sizeof(enum_value_existing));
        existing_enums_set.insert(sv);
    }

    // Find any new enumeration values
    std::vector<ValueType> enum_values_to_add;

    std::optional<std::unordered_set<std::string_view>>
        opt_covered_values = _find_covered_enum_values(
            enum_values_in_write_as_sv, index_schema, index_array);

    if (opt_covered_values.has_value()) {
        const auto& covered_values = opt_covered_values.value();
        const size_t n = enum_values_in_write.size();
        for (size_t i = 0; i < n; i++) {
            const auto& enum_value_in_write = enum_values_in_write[i];
            const auto& sv = enum_values_in_write_as_sv[i];
            if (!existing_enums_set.contains(sv)) {
                if (covered_values.find(sv) != covered_values.end()) {
                    enum_values_to_add.push_back(enum_value_in_write);
                }
            }
        }
    } else {
        const size_t n = enum_values_in_write.size();
        for (size_t i = 0; i < n; i++) {
            const auto& enum_value_in_write = enum_values_in_write[i];
            const auto& sv = enum_values_in_write_as_sv[i];
            if (!existing_enums_set.contains(sv)) {
                enum_values_to_add.push_back(enum_value_in_write);
            }
        }
    }

    // There are two paths to enumeration extension:
    // 1. The user does dataframe.write.
    //    Here, if the on-schema values are a,b,c and the values being
    //    written are c,d,e, then, the expected UX is that values d,e
    //    will be handed to the core enumeration-extend.
    // 2. The user does dataframe.extend_enumeration_values.
    //    Here, if the on-schema values are a,b,c and the values being
    //    extended are c,d,e, then, the expected UX (requested in sc-63930) is:
    //    * By default, that's an error -- they were supposed to just pass d,e.
    //    * If the deduplicate flag was passed, then we allow that, but
    //      we still need to pass core only the d,e values. (It will
    //      throw otherwise).
    if (!deduplicate &&
        enum_values_to_add.size() != enum_values_in_write.size()) {
        throw TileDBSOMAError(fmt::format(
            "[extend_enumeration] one or more values provided are already "
            "present in the enumeration for column '{}', and deduplicate was "
            "not specified",
            column_name));
    }

    // Extend the enumeration in the schema, if there are any values to be
    // added.
    if (enum_values_to_add.size() != 0) {
        // We have new enumeration values; additional processing needed

        // Check if the number of new enumeration will cause an overflow if
        // extended
        auto disk_index_type = schema_->attribute(column_name).type();
        uint64_t max_capacity = _get_max_capacity(disk_index_type);
        auto free_capacity = max_capacity - enum_values_existing.size();
        if (free_capacity < enum_values_to_add.size()) {
            throw TileDBSOMAError(
                "Cannot extend enumeration; reached maximum capacity");
        }

        // Take the existing enumeration values on disk and extend with the new
        // enumeration values
        auto extended_enmr = enmr.extend(enum_values_to_add);
        se.extend_enumeration(extended_enmr);
        return std::tuple{
            true,  // was_extended
            enum_values_in_write,
            extended_enmr};

    } else {
        // The enumeration was not extended
        return std::tuple{
            false,  // was_extended
            enum_values_in_write,
            enmr};
    }
}

std::vector<uint8_t> ManagedQuery::_bool_data_bits_to_bytes(
    ArrowSchema* schema, ArrowArray* array) {
    if (strcmp(schema->format, "b") != 0) {
        throw TileDBSOMAError(fmt::format(
            "_cast_bit_to_uint8 expected column format to be 'b' but saw "
            "{}",
            schema->format));
    }

    const uint8_t* data = reinterpret_cast<const uint8_t*>(array->buffers[1]);
    return *util::bitmap_to_uint8(data, array->length, array->offset);
}

std::optional<std::vector<uint8_t>> ManagedQuery::_cast_validity_buffer(
    ArrowArray* array) {
    const uint8_t* validity = reinterpret_cast<const uint8_t*>(
        array->buffers[0]);
    return util::bitmap_to_uint8(validity, array->length, array->offset);
}

template <>
std::vector<std::string_view> ManagedQuery::_enumeration_values_view(
    Enumeration& enumeration) {
    const void* data;
    uint64_t data_size;

    ctx_->handle_error(tiledb_enumeration_get_data(
        ctx_->ptr().get(), enumeration.ptr().get(), &data, &data_size));

    const void* offsets;
    uint64_t offsets_size;
    ctx_->handle_error(tiledb_enumeration_get_offsets(
        ctx_->ptr().get(), enumeration.ptr().get(), &offsets, &offsets_size));

    std::string_view char_data(static_cast<const char*>(data), data_size);
    const uint64_t* elems = static_cast<const uint64_t*>(offsets);
    size_t count = offsets_size / sizeof(uint64_t);

    std::vector<std::string_view> ret(count);
    for (size_t i = 0; i < count; i++) {
        uint64_t len;
        if (i + 1 < count) {
            len = elems[i + 1] - elems[i];
        } else {
            len = data_size - elems[i];
        }

        ret[i] = char_data.substr(elems[i], len);
    }

    return ret;
}

Enumeration ManagedQuery::get_enumeration(
    std::shared_ptr<Context> ctx,
    std::shared_ptr<Array> arr,
    ArrowSchema* index_schema,
    ArrowSchema* value_schema) {
    return util::get_enumeration(ctx, arr, index_schema, value_schema);
}
};  // namespace tiledbsoma
