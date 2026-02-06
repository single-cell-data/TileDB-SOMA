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
 *   This declares the managed query API.
 */

#include "managed_query.h"

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include "array_buffers.h"
#include "column_buffer.h"
#include "column_buffer_strategies.h"
#include "nanoarrow/nanoarrow.hpp"

#include "../arrow/utils.h"
#include "../arrow/decode.h"

#include "../logging/impl/logger.h"
#include "../logging/logger.h"

namespace tiledbsoma::common {

#pragma region Status implementation

Status::Status()
    : code_(StatusCode::OK) {
}

Status::Status(StatusCode code, std::string_view origin, std::string_view message)
    : code_(code)
    , origin_(origin)
    , message_(message) {
}

#pragma endregion

ManagedQuery::ManagedQuery(
    std::shared_ptr<tiledb::Array> array, std::shared_ptr<tiledb::Context> ctx, std::string_view name)
    : ctx_(ctx)
    , array_(array)
    , name_(name) {
    reset();
}

void ManagedQuery::close() {
    array_->close();
}

void ManagedQuery::reset() {
    query_ = std::make_unique<tiledb::Query>(*ctx_, *array_);
    subarray_ = std::make_unique<tiledb::Subarray>(*ctx_, *array_);

    set_layout(layout_);

    columns_.clear();
    results_complete_ = true;
    total_num_cells_ = 0;
    buffers_.reset();
    query_submitted_ = false;
    retries = 0;
}

void ManagedQuery::set_layout(ResultOrder layout) {
    switch (layout) {
        case ResultOrder::AUTOMATIC:
            if (array_->schema().array_type() == TILEDB_SPARSE)
                query_->set_layout(TILEDB_UNORDERED);
            else
                query_->set_layout(TILEDB_ROW_MAJOR);
            break;
        case ResultOrder::UNORDERED:
            query_->set_layout(TILEDB_UNORDERED);
            break;
        case ResultOrder::GLOBAL:
            query_->set_layout(TILEDB_GLOBAL_ORDER);
            break;
        case ResultOrder::ROWMAJOR:
            query_->set_layout(TILEDB_ROW_MAJOR);
            break;
        case ResultOrder::COLMAJOR:
            query_->set_layout(TILEDB_COL_MAJOR);
            break;
        default:
            throw std::invalid_argument(
                fmt::format("[ManagedQuery][set_layout] invalid ResultOrder({}) passed", static_cast<int>(layout)));
    }

    layout_ = layout;
}

ResultOrder ManagedQuery::result_order() const {
    return layout_;
}

std::shared_ptr<tiledb::Array> ManagedQuery::array() const {
    return array_;
}

std::shared_ptr<ArrayBuffers> ManagedQuery::buffers() const {
    return buffers_;
}

uint64_t ManagedQuery::total_num_cells() const {
    return total_num_cells_;
}

void ManagedQuery::submit_write() {
    query_submitted_ = true;
    setup_write();
    query_->submit();
}

void ManagedQuery::finalize() {
    if (!query_submitted_) {
        throw std::runtime_error(
            "[ManagedQuery] Write query needs to be submitted before "
            "finalizing");
    }
    query_->finalize();
    teardown_write();
}

void ManagedQuery::submit_and_finalize() {
    setup_write();
    query_->submit_and_finalize();
    teardown_write();
}

void ManagedQuery::select_columns(std::span<const std::string> names, bool if_not_empty, bool replace) {
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
        if (!array_->schema().has_attribute(name) && !array_->schema().domain().has_dimension(name)) {
            logging::LOG_WARN(
                fmt::format("[ManagedQuery][select_columns] [{}] Invalid column selected: {}", name_, name));
        } else {
            columns_.push_back(name);
        }
    }
}

void ManagedQuery::reset_columns() {
    columns_.clear();
}

std::vector<std::string> ManagedQuery::column_names() const {
    return columns_;
}

void ManagedQuery::set_condition(const tiledb::QueryCondition& qc) {
    query_->set_condition(qc);
}

void ManagedQuery::set_array_data(ArrowSchema* arrow_schema, ArrowArray* arrow_array) {
    buffers_.reset();

    // Go through all columns in the ArrowTable and cast the values to what is
    // in the ArraySchema on disk
    tiledb::ArraySchemaEvolution se(*ctx_);
    uint64_t ts = array_->open_timestamp_end();
    se.set_timestamp_range(std::make_pair(ts, ts));
    bool evolve_schema = false;
    for (auto i = 0; i < arrow_schema->n_children; ++i) {
        bool enmr_extended = setup_write_arrow_column(arrow_schema->children[i], arrow_array->children[i], se);
        evolve_schema = evolve_schema || enmr_extended;
    }
    if (evolve_schema) {
        se.array_evolve(array_->uri());
    }
}

bool ManagedQuery::setup_write_arrow_column(ArrowSchema* schema, ArrowArray* array, tiledb::ArraySchemaEvolution& se) {
    tiledb_datatype_t disk_type = tiledb_datatype_t::TILEDB_ANY;
    const tiledb::ArraySchema array_schema = array_->schema();

    if (array_schema.has_attribute(schema->name)) {
        std::optional<std::string> enumeration_name = tiledb::AttributeExperimental::get_enumeration_name(
            *ctx_, array_schema.attribute(schema->name));

        if (enumeration_name && schema->dictionary == nullptr) {
            throw std::invalid_argument(
                fmt::format("[ManagedQuery][setup_write_arrow_column] '{}' requires dictionary entry", schema->name));
        }

        // If the attribute is not enumerated, but the provided column is, then we
        // need to use the dictionary values when writing to the array
        if (!enumeration_name && schema->dictionary != nullptr) {
            promote_indexes_to_values(schema, array);
            return false;
        }

        disk_type = array_schema.attribute(schema->name).type();
    } else if (array_schema.domain().has_dimension(schema->name)) {
        disk_type = array_schema.domain().dimension(schema->name).type();
    } else {
        throw std::runtime_error(
            fmt::format(
                "[ManagedQuery][setup_write_arrow_column] TileDB Array is missing a column '{}'", schema->name));
    }

    // In the general cases that do not apply to the two cases above, we need to
    // cast the passed-in column to be what is the type in the schema on disk.
    // Here we identify the passed-in column type (UserType).
    //
    // If _cast_column_aux extended the enumeration then return true. Otherwise,
    // false

    std::unique_ptr<uint8_t[]> validity_buffer = arrow::bitmap_to_bytemap(
        static_cast<const uint8_t*>(array->buffers[0]), array->length, array->offset);

    if (array_schema.has_attribute(schema->name)) {
        const tiledb::Attribute attribute = array_schema.attribute(schema->name);
        std::optional<std::string> enumeration_name = tiledb::AttributeExperimental::get_enumeration_name(
            *ctx_, attribute);

        if (schema->dictionary && enumeration_name) {
            std::optional<std::shared_ptr<tiledb::Enumeration>> extended_enumeration = extend_enumeration(
                schema, array, enumeration_name.value(), se, true);
            remap_enumeration_indices(
                schema,
                array,
                attribute,
                *extended_enumeration.value_or(
                    std::make_shared<tiledb::Enumeration>(
                        tiledb::ArrayExperimental::get_enumeration(*ctx_, *array_, enumeration_name.value().data()))));

            return extended_enumeration.has_value();
        }
    }

    const std::unordered_map<std::string, std::string> metadata = arrow::metadata_string_to_map(schema->metadata);

    if (arrow::to_tiledb_format(schema->format, metadata[]) != disk_type) {
        throw std::runtime_error(
            fmt::format(
                "[ManagedQuery][setup_write_arrow_column] Column '{}'. Expected {}, found {}",
                schema->name,
                tiledb::impl::type_to_str(disk_type),
                tiledb::impl::type_to_str(arrow::to_tiledb_format(schema->format))));
    }

    if (array->n_buffers == 3) {
        if ((strcmp(schema->format, "U") == 0) || (strcmp(schema->format, "Z") == 0)) {
            setup_write_column(
                schema->name,
                array->length,
                static_cast<const std::byte*>(array->buffers[2]) + array->offset,
                static_cast<const uint64_t*>(array->buffers[1]),
                std::move(validity_buffer));
        } else {
            setup_write_column(
                schema->name,
                array->length,
                static_cast<const std::byte*>(array->buffers[2]) + array->offset,
                static_cast<const uint32_t*>(array->buffers[1]),
                std::move(validity_buffer));
        }
    } else {
        if (arrow::to_tiledb_format(schema->format) == TILEDB_BOOL) {
            std::unique_ptr<uint8_t[]> data_bytemap = arrow::bitmap_to_bytemap(
                static_cast<const uint8_t*>(array->buffers[1]), array->length, array->offset);

            setup_write_column(
                schema->name,
                array->length,
                std::unique_ptr<std::byte[]>(reinterpret_cast<std::byte*>(data_bytemap.release())),
                (uint64_t*)nullptr,
                std::move(validity_buffer));
        } else {
            setup_write_column(
                schema->name,
                array->length,
                static_cast<const std::byte*>(array->buffers[1]),
                (uint64_t*)nullptr,
                std::move(validity_buffer));
        }
    }

    return false;
}

std::optional<std::shared_ptr<tiledb::Enumeration>> ManagedQuery::extend_enumeration(
    ArrowSchema* schema,
    ArrowArray* array,
    std::string_view enumeration_name,
    tiledb::ArraySchemaEvolution& se,
    bool deduplicate) {
    tiledb::Enumeration enumeration = tiledb::ArrayExperimental::get_enumeration(
        *ctx_, *array_, enumeration_name.data());

    std::vector<std::string_view> dictionary_values = arrow::dictionary_values_view(
        schema->dictionary, array->dictionary);
    std::unique_ptr<uint8_t[]> validity_buffer = arrow::bitmap_to_bytemap(
        static_cast<const uint8_t*>(array->buffers[0]), array->length, array->offset);

    std::vector<std::string_view> enumeration_values;
    size_t buffer_size = 0;

    auto get_unique_values = [&]<typename IndexType>() {
        // Separate out the values already in the array schema from the values not already in the array schema.
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
        std::span<const IndexType> indices(static_cast<const IndexType*>(array->buffers[1]), array->length);
        std::set<std::string_view> unique_values;

        if (validity_buffer) {
            for (size_t i = 0; i < dictionary_values.size(); ++i) {
                if (validity_buffer[i] && !enumeration.index_of(dictionary_values[i]).has_value() &&
                    !unique_values.contains(dictionary_values[i])) {
                    unique_values.insert(dictionary_values[i]);
                    enumeration_values.push_back(dictionary_values[i]);
                    buffer_size += dictionary_values[i].size();
                }
            }
        } else {
            for (const auto value : dictionary_values) {
                if (!enumeration.index_of(value).has_value() && !unique_values.contains(value)) {
                    unique_values.insert(value);
                    enumeration_values.push_back(value);
                    buffer_size += value.size();
                }
            }
        }
    };

    switch (arrow::to_tiledb_format(schema->format)) {
        case TILEDB_UINT8:
            get_unique_values.template operator()<uint8_t>();
            break;
        case TILEDB_INT8:
            get_unique_values.template operator()<int8_t>();
            break;
        case TILEDB_UINT16:
            get_unique_values.template operator()<uint16_t>();
            break;
        case TILEDB_INT16:
            get_unique_values.template operator()<int16_t>();
            break;
        case TILEDB_UINT32:
            get_unique_values.template operator()<uint32_t>();
            break;
        case TILEDB_INT32:
            get_unique_values.template operator()<int32_t>();
            break;
        case TILEDB_UINT64:
            get_unique_values.template operator()<uint64_t>();
            break;
        case TILEDB_INT64:
            get_unique_values.template operator()<int64_t>();
            break;
        default:
            throw std::runtime_error("");
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
    if (!deduplicate && enumeration_values.size() != dictionary_values.size()) {
        throw std::runtime_error(
            fmt::format(
                "[ManagedQuery][extend_enumeration] one or more values provided are already present in the enumeration "
                "for column "
                "'{}', and deduplicate was not specified",
                schema->name));
    }

    if (enumeration_values.empty()) {
        return std::nullopt;
    }

    // Check if the number of new enumeration will cause an overflow if extended
    size_t free_capacity = get_max_capacity(array_->schema().attribute(schema->name).type()) -
                           enumeration_value_count(*ctx_, enumeration);
    if (free_capacity < enumeration_values.size()) {
        throw std::runtime_error(
            "[ManagedQuery][extend_enumeration] Cannot extend enumeration; reached maximum capacity");
    }

    // Gnerate the buffers to use for extending the enumeration
    std::unique_ptr<std::byte[]> data_buffer = std::make_unique_for_overwrite<std::byte[]>(buffer_size);
    std::optional<std::shared_ptr<tiledb::Enumeration>> extended_enumeration;

    if (enumeration.cell_val_num() == TILEDB_VAR_NUM) {
        std::unique_ptr<uint64_t[]> offsets_buffer = std::make_unique_for_overwrite<uint64_t[]>(
            enumeration_values.size());

        size_t current_offset = 0;
        for (size_t i = 0; i < enumeration_values.size(); ++i) {
            const auto& value = enumeration_values[i];

            std::memcpy(data_buffer.get() + current_offset, value.data(), value.size());
            offsets_buffer[i] = current_offset;
            current_offset += value.size();
        }

        extended_enumeration = std::make_optional(
            std::make_shared<tiledb::Enumeration>(enumeration.extend(
                data_buffer.get(), buffer_size, offsets_buffer.get(), enumeration_values.size() * sizeof(uint64_t))));
    } else {
        size_t stride = tiledb::impl::type_size(enumeration.type()) * enumeration.cell_val_num();
        for (size_t i = 0; i < enumeration_values.size(); ++i) {
            std::memcpy(data_buffer.get() + i * stride, enumeration_values[i].data(), stride);
        }

        extended_enumeration = std::make_optional(
            std::make_shared<tiledb::Enumeration>(enumeration.extend(data_buffer.get(), buffer_size, nullptr, 0)));
    }

    se.extend_enumeration(*extended_enumeration.value());

    return extended_enumeration;
}

std::map<std::string, tiledb::Enumeration> ManagedQuery::attribute_to_enumeration_mapping() const {
    std::map<std::string, tiledb::Enumeration> result;
    const tiledb::ArraySchema schema = array_->schema();
    for (uint32_t i = 0; i < schema.attribute_num(); ++i) {
        tiledb::Attribute attribute = schema.attribute(i);
        std::optional<std::string> enumeration_name = tiledb::AttributeExperimental::get_enumeration_name(
            *ctx_, attribute);

        if (enumeration_name) {
            result.insert(
                {attribute.name(),
                 tiledb::ArrayExperimental::get_enumeration(*ctx_, *array_, enumeration_name.value())});
        }
    }
    return result;
}

void ManagedQuery::remap_enumeration_indices(
    ArrowSchema* schema,
    ArrowArray* array,
    const tiledb::Attribute& attribute,
    const tiledb::Enumeration& enumeration) {
    std::unique_ptr<std::byte[]> data_buffer = std::make_unique_for_overwrite<std::byte[]>(
        array->length * tiledb::impl::type_size(attribute.type()));
    std::unique_ptr<uint8_t[]> validity_buffer = arrow::bitmap_to_bytemap(
        static_cast<const uint8_t*>(array->buffers[0]), array->length, array->offset);

    auto cast_indices = [&]<typename AttributeType, typename IndexType>() {
        std::span<const IndexType> original_indexes(static_cast<const IndexType*>(array->buffers[1]), array->length);
        std::span<AttributeType> shifted_indexes(reinterpret_cast<AttributeType*>(data_buffer.get()), array->length);

        std::vector<std::string_view> values = arrow::dictionary_values_view(schema->dictionary, array->dictionary);
        size_t value_index = 0;
        int value_exist = 0;

        if (validity_buffer && array->null_count != 0) {
            for (size_t i = 0; i < original_indexes.size(); ++i) {
                if (validity_buffer[i] == 1) {
                    tiledb_enumeration_get_value_index(
                        ctx_->ptr().get(),
                        enumeration.ptr().get(),
                        values[original_indexes[i]].data(),
                        values[original_indexes[i]].size(),
                        &value_exist,
                        &value_index);
                    shifted_indexes[i] = static_cast<AttributeType>(value_index);
                } else {
                    shifted_indexes[i] = static_cast<AttributeType>(original_indexes[i]);
                }
            }
        } else {
            for (size_t i = 0; i < original_indexes.size(); ++i) {
                tiledb_enumeration_get_value_index(
                    ctx_->ptr().get(),
                    enumeration.ptr().get(),
                    values[original_indexes[i]].data(),
                    values[original_indexes[i]].size(),
                    &value_exist,
                    &value_index);
                shifted_indexes[i] = static_cast<AttributeType>(value_index);
            }
        }
    };

    auto dispatch_cast = [&]<typename AttributeType>() {
        switch (arrow::to_tiledb_format(schema->format)) {
            case TILEDB_UINT8:
                cast_indices.template operator()<AttributeType, uint8_t>();
                break;
            case TILEDB_INT8:
                cast_indices.template operator()<AttributeType, int8_t>();
                break;
            case TILEDB_UINT16:
                cast_indices.template operator()<AttributeType, uint16_t>();
                break;
            case TILEDB_INT16:
                cast_indices.template operator()<AttributeType, int16_t>();
                break;
            case TILEDB_UINT32:
                cast_indices.template operator()<AttributeType, uint32_t>();
                break;
            case TILEDB_INT32:
                cast_indices.template operator()<AttributeType, int32_t>();
                break;
            case TILEDB_UINT64:
                cast_indices.template operator()<AttributeType, uint64_t>();
                break;
            case TILEDB_INT64:
                cast_indices.template operator()<AttributeType, int64_t>();
                break;
            default:
                throw std::runtime_error(
                    fmt::format(
                        "[ManagedQuery][remap_enumeration_indices][dispatch_cast] Unsupported index type {}",
                        tiledb::impl::type_to_str(arrow::to_tiledb_format(schema->format))));
        }
    };

    switch (attribute.type()) {
        case TILEDB_UINT8:
            dispatch_cast.template operator()<uint8_t>();
            break;
        case TILEDB_INT8:
            dispatch_cast.template operator()<int8_t>();
            break;
        case TILEDB_UINT16:
            dispatch_cast.template operator()<uint16_t>();
            break;
        case TILEDB_INT16:
            dispatch_cast.template operator()<int16_t>();
            break;
        case TILEDB_UINT32:
            dispatch_cast.template operator()<uint32_t>();
            break;
        case TILEDB_INT32:
            dispatch_cast.template operator()<int32_t>();
            break;
        case TILEDB_UINT64:
            dispatch_cast.template operator()<uint64_t>();
            break;
        case TILEDB_INT64:
            dispatch_cast.template operator()<int64_t>();
            break;
        default:
            throw std::runtime_error(
                fmt::format(
                    "[ManagedQuery][remap_enumeration_indices] Unsupported index type {}",
                    tiledb::impl::type_to_str(attribute.type())));
    }

    setup_write_column(
        schema->name, array->length, std::move(data_buffer), (uint64_t*)nullptr, std::move(validity_buffer));
}

void ManagedQuery::promote_indexes_to_values(ArrowSchema* schema, ArrowArray* array) {
    // This is a column with a dictionary. However, the associated TileDB
    // attribute on disk is not enumerated. We will need to map the dictionary
    // indexes to the associated dictionary values and write the values to disk.
    // Here, we identify the passed-in column type

    std::unique_ptr<std::byte[]> data_buffer(nullptr);
    std::unique_ptr<uint64_t[]> offset_buffer(nullptr);
    std::unique_ptr<uint8_t[]> validity_buffer(nullptr);

    switch (arrow::to_tiledb_format(schema->dictionary->format)) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_BLOB:
        case TILEDB_CHAR:
        case TILEDB_GEOM_WKB:
        case TILEDB_GEOM_WKT:
            std::tie(data_buffer, offset_buffer, validity_buffer) = arrow::dictionary_to_values<std::string>(
                schema, array);
            break;
        case TILEDB_BOOL:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<bool>(schema, array);
            break;
        case TILEDB_INT8:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<int8_t>(schema, array);
            break;
        case TILEDB_UINT8:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<uint8_t>(schema, array);
            break;
        case TILEDB_INT16:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<int16_t>(schema, array);
            break;
        case TILEDB_UINT16:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<uint16_t>(schema, array);
            break;
        case TILEDB_INT32:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<int32_t>(schema, array);
            break;
        case TILEDB_UINT32:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<uint32_t>(schema, array);
            break;
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_INT64:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<int64_t>(schema, array);
            break;
        case TILEDB_UINT64:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<uint64_t>(schema, array);
            break;
        case TILEDB_FLOAT32:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<float_t>(schema, array);
            break;
        case TILEDB_FLOAT64:
            std::tie(data_buffer, std::ignore, validity_buffer) = arrow::dictionary_to_values<double_t>(schema, array);
            break;
        default:
            throw std::runtime_error(
                fmt::format(
                    "Saw invalid TileDB value type when attempting to promote indexes to values: {}",
                    tiledb::impl::type_to_str(arrow::to_tiledb_format(schema->dictionary->format))));
    }

    setup_write_column(
        schema->name, array->length, std::move(data_buffer), std::move(offset_buffer), std::move(validity_buffer));
}

void ManagedQuery::setup_write() {
    if (array_->query_type() != TILEDB_WRITE) {
        throw std::runtime_error("[ManagedQuery][setup_write] write requires array to be opened in write mode");
    }

    fill_in_subarrays();

    if (array_->schema().array_type() == TILEDB_DENSE) {
        query_->set_subarray(*subarray_);
    }
}

void ManagedQuery::teardown_write() {
    // Reset
    buffers_.reset();

    // When we evolve the schema, the ArraySchema needs to be updated to the
    // latest version so re-open the Array
    array_->close();
    array_->open(TILEDB_WRITE);

    query_submitted_ = false;
}

void ManagedQuery::setup_read() {
    tiledb::Query::Status status = query_->query_status();

    // If the query is complete nothing to do
    if (status == tiledb::Query::Status::COMPLETE) {
        return;
    }

    fill_in_subarrays();

    // If the query is uninitialized, set the subarray for the query
    if (status == tiledb::Query::Status::UNINITIALIZED) {
        // Set the subarray for range slicing
        query_->set_subarray(*subarray_);
    }

    const tiledb::ArraySchema schema = array_->schema();

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

    if (!buffers_) {
        logging::LOG_TRACE("[ManagedQuery][setup_read] allocate new buffers");
        buffers_ = std::make_shared<ArrayBuffers>(
            columns_, *array_, std::make_unique<MemoryPoolAllocationStrategy>(columns_, *array_));
    }

    for (auto& name : columns_) {
        logging::LOG_DEBUG(
            fmt::format("[ManagedQuery][setup_read] [{}] Attaching buffer for column '{}'", name_, name));
        buffers_->at(name)->attach(*query_);
    }
}

void ManagedQuery::submit_read() {
    query_submitted_ = true;
    query_future_ = std::async(std::launch::async, [&]() {
        logging::LOG_DEBUG("[ManagedQuery] submit thread start");
        try {
            query_->submit();
        } catch (const std::exception& e) {
            return Status(StatusCode::ERROR, "ManagedQuery::submit_read", e.what());
        }
        logging::LOG_DEBUG("[ManagedQuery] submit thread done");
        return Status();
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
        logging::LOG_DEBUG("[ManagedQuery] submit thread start");
        try {
            query_->submit();
        } catch (const std::exception& e) {
            return Status(StatusCode::ERROR, "ManagedQuery::read_next", e.what());
        }
        logging::LOG_DEBUG("[ManagedQuery] submit thread done");
        return Status();
    });

    if (query_future_.valid()) {
        logging::LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Waiting for query", name_));
        query_future_.wait();
        logging::LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Done waiting for query", name_));

        auto retval = query_future_.get();
        if (retval.code() != StatusCode::OK) {
            throw std::runtime_error(fmt::format("[ManagedQuery] [{}] Query FAILED: {}", name_, retval.message()));
        }

    } else {
        throw std::runtime_error(fmt::format("[ManagedQuery] [{}] 'query_future_' invalid", name_));
    }

    tiledb::Query::Status status = query_->query_status();

    if (status == tiledb::Query::Status::FAILED) {
        throw std::runtime_error(fmt::format("[ManagedQuery] [{}] Query FAILED", name_));
    }

    // If the query was ever incomplete, the result buffers contents are not
    // complete.
    if (status == tiledb::Query::Status::INCOMPLETE) {
        results_complete_ = false;
    } else if (status == tiledb::Query::Status::COMPLETE) {
        results_complete_ = true;
    }

    // Update ColumnBuffer size to match query results
    size_t num_cells = 0;
    for (auto& name : buffers_->names()) {
        num_cells = buffers_->at<ReadColumnBuffer>(name)->update_size(*query_);
        logging::LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Buffer {} cells={}", name_, name, num_cells));
    }
    total_num_cells_ += num_cells;

    if (status == tiledb::Query::Status::INCOMPLETE) {
        if (retries > MAX_RETRIES) {
            throw std::runtime_error(fmt::format("[ManagedQuery] [{}] Maximum number or retries reached.", name_));
        }

        tiledb_query_status_details_t status_details;
        tiledb_query_get_status_details(ctx_->ptr().get(), query_->ptr().get(), &status_details);

        switch (status_details.incomplete_reason) {
            case TILEDB_REASON_USER_BUFFER_SIZE:
                // If incomplete because of user buffer size expand buffers otherwise continue submitting the query.
                // Resubmissions without managing to return any data count towards max number or resubmissions.
                if (num_cells == 0) {
                    buffers_->expand_buffers();
                    ++retries;
                }
                break;
            case TILEDB_REASON_MEMORY_BUDGET:
                // If incomplete because of memory budget continue submitting the query.
                // Resubmissions without managing to return any data count towards max number or resubmissions.
                if (num_cells == 0) {
                    ++retries;
                }
                break;
            default:
                ++retries;
        }

        logging::LOG_DEBUG(
            fmt::format(
                "[ManagedQuery] [{}] Query retry; reason {}",
                name_,
                static_cast<int>(status_details.incomplete_reason)));
    }

    return buffers_;
}

void ManagedQuery::fill_in_subarrays() {
    logging::LOG_TRACE("[ManagedQuery] fill_in_subarrays enter");

    // Do this only for dense arrays.
    if (array_->schema().array_type() != TILEDB_DENSE) {
        logging::LOG_TRACE("[ManagedQuery] fill_in_subarrays exit: non-dense");
        return;
    }

    // Don't do this on next-page etc.
    if (query_->query_status() != tiledb::Query::Status::UNINITIALIZED) {
        logging::LOG_TRACE("[ManagedQuery] fill_in_subarrays exit: initialized");
        return;
    }

    auto current_domain = tiledb::ArraySchemaExperimental::current_domain(*ctx_, array_->schema());
    if (current_domain.is_empty()) {
        fill_in_subarrays_domain();
    } else {
        fill_in_subarrays_current_domain(current_domain);
    }
    logging::LOG_TRACE("[ManagedQuery] _fill_in_subarrays exit");
}

void ManagedQuery::fill_in_subarrays_domain() {
    logging::LOG_TRACE("[ManagedQuery][fill_in_subarrays_domain] enter");
    // Dense array must have a subarray set for read. If the array is dense and
    // no ranges have been set, add a range for the array's entire full
    // domain on dimension 0
    if (has_any_subarray_range_set()) {
        return;
    }

    if (array_->query_type() != TILEDB_READ && subarray_->range_num(0) > 0) {
        logging::LOG_TRACE("[ManagedQuery][fill_in_subarrays_domain] range 0 is set");
        return;
    }

    auto [lo, hi] = array_->schema().domain().dimension(0).domain<int64_t>();
    subarray_->add_range(0, lo, hi);
    logging::LOG_TRACE(
        fmt::format(
            "[ManagedQuery][fill_in_subarrays_domain] Add full range to dense subarray dim0 = ({}, {})", lo, hi));

    // Set the subarray for range slicing
    query_->set_subarray(*subarray_);
}

void ManagedQuery::fill_in_subarrays_current_domain(const tiledb::CurrentDomain& current_domain) {
    logging::LOG_TRACE("[ManagedQuery][fill_in_subarrays_current_domain] enter");
    if (current_domain.type() != TILEDB_NDRECTANGLE) {
        throw std::runtime_error("found non-rectangle current-domain type");
    }
    tiledb::NDRectangle ndrect = current_domain.ndrectangle();

    // Loop over dims and apply subarray ranges if not already done by the
    // caller.
    std::vector<tiledb::Dimension> dimensions = array_->schema().domain().dimensions();
    for (size_t i = 0; i < dimensions.size(); ++i) {
        std::string dim_name = dimensions[i].name();
        tiledb_datatype_t dim_type = dimensions[i].type();

        if (subarray_range_set_[dim_name]) {
            logging::LOG_TRACE(fmt::format("[ManagedQuery][fill_in_subarrays_current_domain] continue {}", dim_name));
            continue;
        }

        // Dense arrays are (as of this writing in 1.15.0) all DenseNDArray.
        // Per the spec DenseNDArray must only have dims named
        // soma_dim_{i} with i=0,1,2,...,n-1, of type int64.
        if (dim_name.rfind("soma_dim_", 0) != 0) {
            throw std::runtime_error(fmt::format("found dense array with unexpected dim name {}", dim_name));
        }
        if (dim_type != TILEDB_INT64) {
            throw std::runtime_error(
                fmt::format(
                    "expected dense arrays to have int64 dims; got {} for {}",
                    tiledb::impl::to_str(dim_type),
                    dim_name));
        }

        auto [cd_lo, cd_hi] = ndrect.range<int64_t>(dim_name);
        int64_t lo = cd_lo;
        int64_t hi = cd_hi;

        logging::LOG_TRACE(
            fmt::format(
                "[ManagedQuery] fill_in_subarrays_current_domain dim "
                "name {} current domain ({}, {})",
                dim_name,
                cd_lo,
                cd_hi));

        if (array_->query_type() == TILEDB_READ) {
            auto [dom_lo, dom_hi] = dimensions[i].domain<int64_t>();
            logging::LOG_TRACE(
                fmt::format(
                    "[ManagedQuery] fill_in_subarrays_current_domain "
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

        logging::LOG_TRACE(
            fmt::format(
                "[ManagedQuery] fill_in_subarrays_current_domain dim "
                "name {} select ({}, {})",
                dim_name,
                lo,
                hi));

        select_ranges<int64_t>(dim_name, {{std::make_pair(lo, hi)}});
    }
}

bool ManagedQuery::has_any_empty_range() const {
    for (auto subdim : subarray_range_empty_) {
        if (subdim.second == true) {
            return true;
        }
    }
    return false;
}

bool ManagedQuery::has_any_subarray_range_set() const {
    for (auto subdim : subarray_range_set_) {
        if (subdim.second == true) {
            return true;
        }
    }
    return false;
}

bool ManagedQuery::is_empty_query() const {
    return has_any_empty_range() && has_any_subarray_range_set();
}

bool ManagedQuery::is_complete(bool query_status_only) const {
    return query_->query_status() == tiledb::Query::Status::COMPLETE || (!query_status_only && is_empty_query());
}

bool ManagedQuery::results_complete() const {
    return is_complete() && results_complete_;
}
}  // namespace tiledbsoma::common
