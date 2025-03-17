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
#include <format>
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
            throw std::invalid_argument(std::format(
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
            LOG_WARN(std::format(
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
        LOG_DEBUG(std::format(
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
        LOG_DEBUG(std::format("[ManagedQuery] [{}] Waiting for query", name_));
        query_future_.wait();
        LOG_DEBUG(
            std::format("[ManagedQuery] [{}] Done waiting for query", name_));

        auto retval = query_future_.get();
        if (!retval.succeeded()) {
            throw TileDBSOMAError(std::format(
                "[ManagedQuery] [{}] Query FAILED: {}",
                name_,
                retval.message()));
        }

    } else {
        throw TileDBSOMAError(
            std::format("[ManagedQuery] [{}] 'query_future_' invalid", name_));
    }

    auto status = query_->query_status();

    if (status == Query::Status::FAILED) {
        throw TileDBSOMAError(
            std::format("[ManagedQuery] [{}] Query FAILED", name_));
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
        LOG_DEBUG(std::format(
            "[ManagedQuery] [{}] Buffer {} cells={}", name_, name, num_cells));
    }
    total_num_cells_ += num_cells;

    // TODO: retry the query with larger buffers
    if (status == Query::Status::INCOMPLETE && !num_cells) {
        throw TileDBSOMAError(
            std::format("[ManagedQuery] [{}] Buffers are too small.", name_));
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
    LOG_TRACE(std::format(
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
            LOG_TRACE(std::format(
                "[ManagedQuery] _fill_in_subarrays continue {}", dim_name));
            continue;
        }

        // Dense arrays are (as of this writing in 1.15.0) all DenseNDArray.
        // Per the spec DenseNDArray must only have dims named
        // soma_dim_{i} with i=0,1,2,...,n-1, of type int64.
        if (dim_name.rfind("soma_dim_", 0) != 0) {
            throw TileDBSOMAError(std::format(
                "found dense array with unexpected dim name {}", dim_name));
        }
        if (dim.type() != TILEDB_INT64) {
            throw TileDBSOMAError(std::format(
                "expected dense arrays to have int64 dims; got {} for {}",
                tiledb::impl::to_str(dim.type()),
                dim_name));
        }

        std::array<int64_t, 2> lo_hi_arr = ndrect.range<int64_t>(dim_name);
        int64_t cd_lo = lo_hi_arr[0];
        int64_t cd_hi = lo_hi_arr[1];
        int64_t lo = cd_lo;
        int64_t hi = cd_hi;

        LOG_TRACE(std::format(
            "[ManagedQuery] _fill_in_subarrays_if_dense_with_new_shape dim "
            "name {} current domain ({}, {})",
            dim_name,
            cd_lo,
            cd_hi));

        if (is_read) {
            auto [dom_lo, dom_hi] = schema.domain()
                                        .dimension(0)
                                        .domain<int64_t>();
            LOG_TRACE(std::format(
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

        LOG_TRACE(std::format(
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
        LOG_DEBUG(std::format("[ManagedQuery] [{}] Waiting for query", name_));
        query_future_.wait();
        LOG_DEBUG(
            std::format("[ManagedQuery] [{}] Done waiting for query", name_));

        auto retval = query_future_.get();
        if (!retval.succeeded()) {
            throw TileDBSOMAError(std::format(
                "[ManagedQuery] [{}] Query FAILED: {}",
                name_,
                retval.message()));
        }

    } else {
        throw TileDBSOMAError(
            std::format("[ManagedQuery] [{}] 'query_future_' invalid", name_));
    }

    auto status = query_->query_status();

    if (status == Query::Status::FAILED) {
        throw TileDBSOMAError(
            std::format("[ManagedQuery] [{}] Query FAILED", name_));
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
        LOG_DEBUG(std::format(
            "[ManagedQuery] [{}] Buffer {} cells={}", name_, name, num_cells));
    }
    total_num_cells_ += num_cells;

    // TODO: retry the query with larger buffers
    if (status == Query::Status::INCOMPLETE && !num_cells) {
        throw TileDBSOMAError(
            std::format("[ManagedQuery] [{}] Buffers are too small.", name_));
    }

    return buffers_;
}

void ManagedQuery::check_column_name(const std::string& name) {
    if (!buffers_->contains(name)) {
        throw TileDBSOMAError(std::format(
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
            throw TileDBSOMAError(std::format(
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
            throw TileDBSOMAError(std::format(
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
    std::vector<uint8_t> values = _cast_bool_data(value_schema, value_array);
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
        throw TileDBSOMAError(std::format(
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

    auto casted = _cast_bool_data(schema, array);

    setup_write_column(
        schema->name,
        array->length,
        (const void*)casted.data(),
        (uint64_t*)nullptr,
        _cast_validity_buffer(array));
    return false;
}

bool ManagedQuery::_extend_enumeration(
    ArrowSchema* value_schema,
    ArrowArray* value_array,
    ArrowSchema* index_schema,
    ArrowArray* index_array,
    ArraySchemaEvolution se) {
    // For columns with dictionaries, we need to identify the data type of the
    // enumeration to extend any new enumeration values

    // As of 1.16, and enumeration names are in the format of
    // {index name}_{value dtype}. Prior to 1.16, enum labels are the same as
    // the index name. If the new format doesn't work then fall back to the old
    // format
    Enumeration enmr = [&]() {
        try {
            return ArrayExperimental::get_enumeration(
                *ctx_,
                *array_,
                util::get_enmr_label(index_schema, value_schema));
        } catch (const std::exception& e) {
            return ArrayExperimental::get_enumeration(
                *ctx_, *array_, index_schema->name);
        }
    }();
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
            throw TileDBSOMAError(std::format(
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
        auto casted = _cast_bool_data(value_schema, value_array);
        enums_in_write.assign(casted.data(), casted.data() + num_elems);
    } else {
        // General case
        ValueType* data;
        if (value_array->n_buffers == 3) {
            data = (ValueType*)value_array->buffers[2] + value_array->offset;
        } else {
            data = (ValueType*)value_array->buffers[1] + value_array->offset;
        }
        enums_in_write.assign((ValueType*)data, (ValueType*)data + num_elems);
    }

    // Get all the enumeration values in the on-disk TileDB attribute

    // As of 1.16, and enumeration names are in the format of
    // {index name}_{value dtype}. Prior to 1.16, enum labels are the same as
    // the index name. If the new format doesn't work then fall back to the old
    // format
    Enumeration enmr = [&]() {
        try {
            return ArrayExperimental::get_enumeration(
                *ctx_,
                *array_,
                util::get_enmr_label(index_schema, value_schema));
        } catch (const std::exception& e) {
            return ArrayExperimental::get_enumeration(
                *ctx_, *array_, index_schema->name);
        }
    }();
    std::vector<ValueType> enums_existing = enmr.as_vector<ValueType>();

    // Find any new enumeration values
    std::vector<ValueType> extend_values;
    for (auto enum_val : enums_in_write) {
        if (std::find(enums_existing.begin(), enums_existing.end(), enum_val) ==
            enums_existing.end()) {
            extend_values.push_back(enum_val);
        }
    }

    std::string column_name = index_schema->name;
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

    std::string_view data(
        static_cast<const char*>(value_array->buffers[2]), offsets_v.back());

    std::vector<std::string_view> enums_in_write;
    for (size_t i = 0; i < num_elems; ++i) {
        auto beg = offsets_v[i];
        auto sz = offsets_v[i + 1] - beg;
        enums_in_write.push_back(data.substr(beg, sz));
    }

    // As of 1.16, and enumeration names are in the format of
    // {index name}_{value dtype}. Prior to 1.16, enum labels are the same as
    // the index name. If the new format doesn't work then fall back to the old
    // format
    Enumeration enmr = [&]() {
        try {
            return ArrayExperimental::get_enumeration(
                *ctx_,
                *array_,
                util::get_enmr_label(index_schema, value_schema));
        } catch (const std::exception& e) {
            return ArrayExperimental::get_enumeration(
                *ctx_, *array_, index_schema->name);
        }
    }();
    std::vector<std::string_view> extend_values;
    size_t total_size = 0;

    auto enums_existing = _enumeration_values_view<std::string_view>(enmr);
    std::unordered_set<std::string_view> existing_enums_set;
    for (const auto& existing_enum_val : enums_existing) {
        existing_enums_set.insert(existing_enum_val);
    }

    for (const auto& enum_val : enums_in_write) {
        if (!existing_enums_set.contains(enum_val)) {
            extend_values.push_back(enum_val);
            total_size += enum_val.size();
        }
    }

    std::string column_name = index_schema->name;
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

        std::vector<uint8_t> extend_data(total_size);
        std::vector<uint64_t> extend_offsets;
        extend_offsets.reserve(extend_values.size());
        size_t curr_offset = 0;
        for (const auto& enum_val : extend_values) {
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
            extend_values.size() * sizeof(uint64_t));
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

std::vector<uint8_t> ManagedQuery::_cast_bool_data(
    ArrowSchema* schema, ArrowArray* array) {
    if (strcmp(schema->format, "b") != 0) {
        throw TileDBSOMAError(std::format(
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
};  // namespace tiledbsoma
