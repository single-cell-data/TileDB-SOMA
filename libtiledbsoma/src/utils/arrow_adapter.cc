/**
 * @file   arrow_adapter.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the ArrowAdapter class.
 */

#include <ranges>

#include "../soma/column_buffer.h"
#include "arrow_adapter.h"
#include "logger.h"
#include "util.h"

#include "../soma/soma_attribute.h"
#include "../soma/soma_coordinates.h"
#include "../soma/soma_dimension.h"
#include "../soma/soma_geometry_column.h"

namespace tiledbsoma {

using namespace tiledb;

/**************************
 * Internal helper methods
 **************************/

std::string_view to_arrow_readable(std::string_view arrow_dtype) {
    std::map<std::string_view, std::string_view> _to_arrow_readable = {
        {"n", "null"},
        {"b", "boolean"},
        {"c", "int8"},
        {"C", "uint8"},
        {"s", "int16"},
        {"S", "uint16"},
        {"i", "int32"},
        {"I", "uint32"},
        {"l", "int64"},
        {"L", "uint64"},
        {"e", "float16"},
        {"f", "float32"},
        {"g", "float64"},
        {"z", "binary"},
        {"Z", "large binary"},
        {"vz", "binary view"},
        {"u", "utf-8 string"},
        {"U", "large utf-8 string"},
        {"vu", "utf-8 view"},
        {"tdD", "date32 [days]"},
        {"tdm", "date64 [milliseconds]"},
        {"tts", "time32 [seconds]"},
        {"ttm", "time32 [milliseconds]"},
        {"ttu", "time64 [microseconds]"},
        {"ttn", "time64 [nanoseconds]"},
        {"tDs", "duration [seconds]"},
        {"tDm", "duration [milliseconds]"},
        {"tDu", "duration [microseconds]"},
        {"tDn", "duration [nanoseconds]"},
        {"tiM", "interval [months]"},
        {"tiD", "interval [days, time]"},
        {"tin", "interval [month, day, nanoseconds]"},
        {"+l", "list"},
        {"+L", "large list"},
        {"+vl", "list-view"},
        {"+vL", "large list-view"},
        {"+s", "struct"},
        {"+m", "map"},
        {"+r", "run-end encoded"}};

    auto it = _to_arrow_readable.find(arrow_dtype);
    return it != _to_arrow_readable.end() ? it->second :
                                            "unknown Arrow type [see "
                                            "https://arrow.apache.org/docs/format/"
                                            "CDataInterface.html#data-type-description-format-strings]";
}

enum ArrowType to_nanoarrow_type(std::string_view arrow_dtype) {
    std::map<std::string_view, enum ArrowType> _to_nanoarrow_type_map = {
        {"i", NANOARROW_TYPE_INT32},        {"c", NANOARROW_TYPE_INT8},         {"C", NANOARROW_TYPE_UINT8},
        {"s", NANOARROW_TYPE_INT16},        {"S", NANOARROW_TYPE_UINT16},       {"I", NANOARROW_TYPE_UINT32},
        {"l", NANOARROW_TYPE_INT64},        {"L", NANOARROW_TYPE_UINT64},       {"f", NANOARROW_TYPE_FLOAT},
        {"g", NANOARROW_TYPE_DOUBLE},       {"u", NANOARROW_TYPE_STRING},       {"U", NANOARROW_TYPE_LARGE_STRING},
        {"b", NANOARROW_TYPE_BOOL},         {"tss:", NANOARROW_TYPE_TIMESTAMP}, {"tsm:", NANOARROW_TYPE_TIMESTAMP},
        {"tsn:", NANOARROW_TYPE_TIMESTAMP}, {"tsu:", NANOARROW_TYPE_TIMESTAMP}, {"tdD", NANOARROW_TYPE_TIMESTAMP},
        {"z", NANOARROW_TYPE_BINARY},       {"Z", NANOARROW_TYPE_LARGE_BINARY},
    };

    try {
        return _to_nanoarrow_type_map.at(arrow_dtype);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(
            fmt::format("ArrowAdapter: Unsupported Arrow type: {} ({})", arrow_dtype, to_arrow_readable(arrow_dtype)));
    }
}

std::pair<enum ArrowType, enum ArrowTimeUnit> to_nanoarrow_time(std::string_view arrow_dtype) {
    std::map<std::string_view, std::pair<enum ArrowType, enum ArrowTimeUnit>> _to_nanoarrow_time = {
        {"tss:", {NANOARROW_TYPE_TIMESTAMP, NANOARROW_TIME_UNIT_SECOND}},
        {"tsm:", {NANOARROW_TYPE_TIMESTAMP, NANOARROW_TIME_UNIT_MILLI}},
        {"tsu:", {NANOARROW_TYPE_TIMESTAMP, NANOARROW_TIME_UNIT_MICRO}},
        {"tsn:", {NANOARROW_TYPE_TIMESTAMP, NANOARROW_TIME_UNIT_NANO}},
    };

    try {
        return _to_nanoarrow_time.at(arrow_dtype);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(
            fmt::format("ArrowAdapter: Unsupported Arrow type: {} ({})", arrow_dtype, to_arrow_readable(arrow_dtype)));
    }
}

ArrowBuffer::ArrowBuffer(ReadColumnBuffer* buffer, bool large_offsets) {
    if (buffer->is_var()) {
        size_t data_byte_size = buffer->offsets()[buffer->size()];
        data_ = std::make_unique_for_overwrite<std::byte[]>(data_byte_size);

        if (large_offsets) {
            size_t offset_byte_size = (buffer->size() + 1) * sizeof(int64_t);
            large_offsets_ = std::make_unique_for_overwrite<int64_t[]>(buffer->size() + 1);
            std::memcpy(large_offsets_.get(), buffer->offsets().data(), offset_byte_size);
        } else {
            small_offsets_ = std::make_unique_for_overwrite<int32_t[]>(buffer->is_var() + 1);
            auto offsets = buffer->offsets();
            for (size_t i = 0; i < offsets.size(); ++i) {
                small_offsets_[i] = static_cast<int32_t>(offsets[i]);
            }
        }

        std::memcpy(data_.get(), buffer->data().data(), data_byte_size);
    } else {
        if (buffer->type() == TILEDB_BOOL) {
            size_t data_byte_size = (buffer->size() + 7) / 8;

            data_ = std::make_unique_for_overwrite<std::byte[]>(data_byte_size);
            buffer->data_to_bitmap();

            std::memcpy(data_.get(), buffer->data().data(), data_byte_size);
        } else {
            size_t data_byte_size = buffer->size() * tiledb::impl::type_size(buffer->type());

            data_ = std::make_unique_for_overwrite<std::byte[]>(data_byte_size);

            std::memcpy(data_.get(), buffer->data().data(), data_byte_size);
        }
    }

    if (buffer->is_nullable()) {
        buffer->validity_to_bitmap();
        auto bitmap_size = (buffer->size() + 7) / 8;

        validity_ = std::make_unique_for_overwrite<std::byte[]>(bitmap_size);
        std::memcpy(validity_.get(), buffer->validity().data(), bitmap_size);
    }

    length = buffer->size();
    name = buffer->name();
}

ArrowBuffer::ArrowBuffer(const Enumeration& enumeration, bool large_offsets) {
    Context ctx = enumeration.context();

    const void* data;
    uint64_t data_size;
    ctx.handle_error(tiledb_enumeration_get_data(ctx.ptr().get(), enumeration.ptr().get(), &data, &data_size));

    data_ = std::make_unique_for_overwrite<std::byte[]>(data_size);
    std::memcpy(data_.get(), data, data_size);

    switch (enumeration.type()) {
        case TILEDB_CHAR:
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_BLOB:
        case TILEDB_GEOM_WKT:
        case TILEDB_GEOM_WKB: {
            const void* offsets;
            uint64_t offsets_size;
            ctx.handle_error(
                tiledb_enumeration_get_offsets(ctx.ptr().get(), enumeration.ptr().get(), &offsets, &offsets_size));
            size_t count = offsets_size / sizeof(uint64_t);

            if (large_offsets) {
                large_offsets_ = std::make_unique_for_overwrite<int64_t[]>(count + 1);
                std::memcpy(large_offsets_.get(), offsets, offsets_size);
                large_offsets_[count] = data_size;
            } else {
                small_offsets_ = std::make_unique_for_overwrite<int32_t[]>(count + 1);
                std::span<const uint64_t> offsets_v(static_cast<const uint64_t*>(offsets), count);
                for (size_t i = 0; i < count; ++i) {
                    small_offsets_[i] = static_cast<int32_t>(offsets_v[i]);
                }
                small_offsets_[count] = static_cast<int32_t>(data_size);
            }

            length = count;
        } break;
        case TILEDB_BOOL: {
            data_ = std::make_unique_for_overwrite<std::byte[]>(1);
            std::span<const bool> data_v(static_cast<const bool*>(data), data_size);
            size_t count = data_size / sizeof(bool);

            // Represent the Boolean vector with, at most, the last two
            // bits. In Arrow, Boolean values are LSB packed
            uint8_t packed_data = 0;
            for (size_t i = 0; i < count; ++i)
                packed_data |= (data_v[i] << i);

            std::memcpy(data_.get(), &packed_data, 1);
            length = count;
        } break;
        case TILEDB_INT8:
            length = data_size / sizeof(int8_t);
            break;
        case TILEDB_UINT8:
            length = data_size / sizeof(uint8_t);
            break;
        case TILEDB_INT16:
            length = data_size / sizeof(int16_t);
            break;
        case TILEDB_UINT16:
            length = data_size / sizeof(uint16_t);
            break;
        case TILEDB_INT32:
            length = data_size / sizeof(int32_t);
            break;
        case TILEDB_UINT32:
            length = data_size / sizeof(uint32_t);
            break;
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_INT64:
            length = data_size / sizeof(int64_t);
            break;
        case TILEDB_UINT64:
            length = data_size / sizeof(uint64_t);
            break;
        case TILEDB_FLOAT32:
            length = data_size / sizeof(float_t);
            break;
        case TILEDB_FLOAT64:
            length = data_size / sizeof(double_t);
            break;
        default:
            throw TileDBSOMAError(
                fmt::format(
                    "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                    tiledb::impl::type_to_str(enumeration.type())));
    }

    name = enumeration.name();
}

/**************************
 * External API
 **************************/

void ArrowAdapter::release_schema(struct ArrowSchema* schema) {
    std::string name_for_log(schema->name == nullptr ? "anonymous" : schema->name);
    if (schema->name != nullptr)
        LOG_DEBUG(fmt::format("[ArrowAdapter] release_schema start for {}", schema->name));

    if (schema->name != nullptr) {
        LOG_TRACE(fmt::format("[ArrowAdapter] release_schema schema->name {}", schema->name));
        free((void*)schema->name);
        schema->name = nullptr;
    }
    if (schema->format != nullptr) {
        LOG_TRACE(fmt::format("[ArrowAdapter] release_schema name {} schema->format {}", name_for_log, schema->format));
        free((void*)schema->format);
        schema->format = nullptr;
    }
    if (schema->metadata != nullptr) {
        LOG_TRACE(fmt::format("[ArrowAdapter] release_schema name {} schema->metadata", name_for_log));
        free((void*)schema->metadata);
        schema->metadata = nullptr;
    }

    if (schema->children != nullptr) {
        LOG_TRACE(
            fmt::format(
                "[ArrowAdapter] release_schema name {} n_children {} begin "
                "recurse ",
                name_for_log,
                schema->n_children));

        for (auto i = 0; i < schema->n_children; i++) {
            if (schema->children[i] != nullptr) {
                if (schema->children[i]->release != nullptr) {
                    LOG_TRACE(
                        fmt::format(
                            "[ArrowAdapter] release_schema name {} schema->child "
                            "{} "
                            "release",
                            name_for_log,
                            i));
                    schema->children[i]->release(schema->children[i]);
                }
                LOG_TRACE(
                    fmt::format(
                        "[ArrowAdapter] release_schema name {} schema->child {} "
                        "free",
                        name_for_log,
                        i));
                free(schema->children[i]);
                schema->children[i] = nullptr;
            }
        }

        LOG_TRACE(
            fmt::format(
                "[ArrowAdapter] release_schema name {} n_children {} end recurse ", name_for_log, schema->n_children));

        free(schema->children);
        schema->children = nullptr;
    }

    if (schema->dictionary != nullptr) {
        if (schema->dictionary->release != nullptr) {
            LOG_TRACE(fmt::format("[ArrowAdapter] release_schema name {} schema->dict release", name_for_log));
            schema->dictionary->release(schema->dictionary);
        }
        LOG_TRACE(fmt::format("[ArrowAdapter] release_schema name {} schema->dict free", name_for_log));
        free(schema->dictionary);
        schema->dictionary = nullptr;
    }

    schema->release = nullptr;
    LOG_TRACE(fmt::format("[ArrowAdapter] release_schema name {} done", name_for_log));
}

void ArrowAdapter::release_array(struct ArrowArray* array) {
    auto arrow_buffer = static_cast<PrivateArrowBuffer*>(array->private_data);
    if (arrow_buffer != nullptr) {
        LOG_TRACE(
            fmt::format(
                "[ArrowAdapter] release_array {} use_count={}",
                arrow_buffer->buffer_->name,
                arrow_buffer->buffer_.use_count()));

        // Delete the ArrowBuffer, which was allocated with new.
        // If the ArrowBuffer.buffer_ shared_ptr is the last reference to the
        // underlying ColumnBuffer, the ColumnBuffer will be deleted.
        delete arrow_buffer;
    } else {
        // ArrowArray that allocate buffers directly to ``array->buffers`` will
        // leak memory. We detect these arrays when the ``array->private_data``
        // is null
        for (int64_t i = 0; i < array->n_buffers; ++i) {
            if (array->buffers[i] != nullptr) {
                free((void*)array->buffers[i]);
                array->buffers[i] = nullptr;
            }
        }
    }

    if (array->buffers != nullptr) {
        free(array->buffers);
        array->buffers = nullptr;
    }

    if (array->children != nullptr) {
        for (auto i = 0; i < array->n_children; i++) {
            if (array->children[i] != nullptr) {
                if (array->children[i]->release != nullptr) {
                    LOG_TRACE(fmt::format("[ArrowAdapter] release_schema array->child {} release", i));

                    array->children[i]->release(array->children[i]);
                }
                LOG_TRACE(fmt::format("[ArrowAdapter] release_schema array->child {} free", i));
                free(array->children[i]);
                array->children[i] = nullptr;
            }
        }
        LOG_TRACE("[ArrowAdapter] release_array array->children");
        free(array->children);
        array->children = nullptr;
    }

    if (array->dictionary != nullptr) {
        LOG_TRACE("[ArrowAdapter] release_array array->dict release");

        array->dictionary->release(array->dictionary);
        free(array->dictionary);
        array->dictionary = nullptr;
    }

    array->release = nullptr;
    LOG_TRACE(fmt::format("[ArrowAdapter] release_array done"));
}

ArrowSchema* ArrowAdapter::arrow_schema_from_tiledb_dimension(const Dimension& dimension) {
    // Accessing dimension attributes may throw.
    // To avoid leaking memory we need to access the before allocating any
    // memory so we are able to free any allocated memory of the caller without
    // any leaking from this function
    auto format = ArrowAdapter::to_arrow_format(dimension.type());
    auto name = dimension.name();

    ArrowSchema* arrow_schema = (ArrowSchema*)malloc(sizeof(ArrowSchema));
    arrow_schema->format = strdup(format.data());
    arrow_schema->name = strdup(name.c_str());
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    arrow_schema->n_children = 0;
    arrow_schema->children = nullptr;
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;
    LOG_TRACE(
        fmt::format(
            "[ArrowAdapter] arrow_schema_from_tiledb_dimension format {} "
            "name {}",
            arrow_schema->format,
            arrow_schema->name));

    return arrow_schema;
}

ArrowSchema* ArrowAdapter::arrow_schema_from_tiledb_attribute(
    const Attribute& attribute, const Context& ctx, const Array& tiledb_array, bool downcast_dict_of_large_var) {
    // Accessing dimension attributes may throw.
    // To avoid leaking memory we need to access the before allocating any
    // memory so we are able to free any allocated memory of the caller without
    // any leaking from this function
    auto format = ArrowAdapter::to_arrow_format(attribute.type());
    auto name = attribute.name();
    auto enmr_name = AttributeExperimental::get_enumeration_name(ctx, attribute);
    auto enmr = enmr_name ? std::make_optional(ArrayExperimental::get_enumeration(ctx, tiledb_array, *enmr_name)) :
                            std::nullopt;

    ArrowSchema* arrow_schema = (ArrowSchema*)malloc(sizeof(ArrowSchema));
    arrow_schema->format = strdup(format.data());
    arrow_schema->name = strdup(name.c_str());
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    if (attribute.nullable() && attribute.name() != SOMA_GEOMETRY_COLUMN_NAME) {
        arrow_schema->flags |= ARROW_FLAG_NULLABLE;
    } else {
        arrow_schema->flags &= ~ARROW_FLAG_NULLABLE;
    }
    arrow_schema->n_children = 0;
    arrow_schema->children = nullptr;
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;

    if (attribute.type() == TILEDB_GEOM_WKB) {
        nanoarrow::UniqueBuffer metadata_buffer;
        ArrowMetadataBuilderInit(metadata_buffer.get(), nullptr);
        ArrowMetadataBuilderAppend(metadata_buffer.get(), ArrowCharView("dtype"), ArrowCharView("WKB"));
        ArrowSchemaSetMetadata(arrow_schema, reinterpret_cast<char*>(metadata_buffer->data));
    }

    LOG_TRACE(
        fmt::format(
            "[ArrowAdapter] arrow_schema_from_tiledb_array format {} "
            "name {}",
            arrow_schema->format,
            arrow_schema->name));

    if (enmr_name.has_value()) {
        auto dict = (ArrowSchema*)malloc(sizeof(ArrowSchema));
        if (downcast_dict_of_large_var) {
            if (enmr->type() == TILEDB_STRING_ASCII || enmr->type() == TILEDB_CHAR) {
                dict->format = strdup("z");
            } else {
                dict->format = strdup(ArrowAdapter::to_arrow_format(enmr->type(), false).data());
            }
        } else {
            dict->format = strdup(ArrowAdapter::to_arrow_format(enmr->type()).data());
        }
        dict->name = strdup(enmr->name().c_str());
        dict->metadata = nullptr;
        if (enmr->ordered()) {
            arrow_schema->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
        } else {
            arrow_schema->flags &= ~ARROW_FLAG_DICTIONARY_ORDERED;
        }
        dict->n_children = 0;
        dict->children = nullptr;
        dict->dictionary = nullptr;
        dict->release = &ArrowAdapter::release_schema;
        dict->private_data = nullptr;
        arrow_schema->dictionary = dict;
    }

    return arrow_schema;
}

std::tuple<ArraySchema, nlohmann::json> ArrowAdapter::tiledb_schema_from_arrow_schema(
    std::shared_ptr<Context> ctx,
    const managed_unique_ptr<ArrowSchema>& arrow_schema,
    const ArrowTable& index_column_info,
    const std::optional<SOMACoordinateSpace>& coordinate_space,
    std::string soma_type,
    bool is_sparse,
    PlatformConfig platform_config,
    std::optional<std::pair<int64_t, int64_t>> timestamp_range) {
    auto schema = utils::create_base_tiledb_schema(ctx, platform_config, is_sparse, timestamp_range);

    auto& index_column_array = index_column_info.first;
    auto& index_column_schema = index_column_info.second;

    std::vector<std::shared_ptr<SOMAColumn>> columns;

    for (int64_t sch_idx = 0; sch_idx < arrow_schema->n_children; ++sch_idx) {
        auto child = arrow_schema->children[sch_idx];
        std::string_view type_metadata;

        if (ArrowMetadataHasKey(child->metadata, ArrowCharView(ARROW_DATATYPE_METADATA_KEY.c_str()))) {
            ArrowStringView out;
            NANOARROW_THROW_NOT_OK(
                ArrowMetadataGetValue(child->metadata, ArrowCharView(ARROW_DATATYPE_METADATA_KEY.c_str()), &out));

            type_metadata = std::string_view(out.data, out.size_bytes);
        }

        LOG_DEBUG(fmt::format("[ArrowAdapter] schema pass for child {} name '{}'", sch_idx, std::string(child->name)));

        bool isattr = true;

        for (int64_t i = 0; i < index_column_schema->n_children; ++i) {
            if (strcmp(child->name, index_column_schema->children[i]->name) == 0) {
                if (strcmp(child->name, SOMA_GEOMETRY_COLUMN_NAME.c_str()) == 0) {
                    columns.push_back(
                        SOMAGeometryColumn::create(
                            ctx,
                            child,
                            index_column_schema->children[i],
                            index_column_array->children[i],
                            coordinate_space.value(),
                            soma_type,
                            type_metadata,
                            platform_config));
                } else {
                    columns.push_back(
                        SOMADimension::create(
                            ctx,
                            index_column_schema->children[i],
                            index_column_array->children[i],
                            soma_type,
                            type_metadata,
                            platform_config));
                }
                isattr = false;
                LOG_DEBUG(fmt::format("[ArrowAdapter] adding dimension {}", child->name));
                break;
            }
        }

        if (isattr) {
            columns.push_back(SOMAAttribute::create(ctx, child, type_metadata, platform_config));
            LOG_DEBUG(fmt::format("[ArrowAdapter] adding attribute {}", child->name));
        }
    }

    LOG_DEBUG(fmt::format("[ArrowAdapter] Additional schema metadata"));
    nlohmann::json soma_schema_extension;
    soma_schema_extension[TILEDB_SOMA_SCHEMA_COL_KEY] = nlohmann::json::array();
    soma_schema_extension["version"] = TILEDB_SOMA_SCHEMA_VERSION;

    Domain domain(*ctx);
    // Unit tests expect dimension order should match the index column schema
    // and NOT the Arrow schema
    // We generate the additional schema metadata here to ensure that the
    // serialized column order matches the expected schema order
    for (int64_t i = 0; i < index_column_schema->n_children; ++i) {
        LOG_DEBUG(fmt::format("[ArrowAdapter] child {}", i));
        const auto column = util::find_column_by_name(columns, index_column_schema->children[i]->name);

        column->serialize(soma_schema_extension[TILEDB_SOMA_SCHEMA_COL_KEY]);

        if (column->tiledb_dimensions().has_value()) {
            // Intermediate variable required to avoid lifetime issues
            auto dimensions = column->tiledb_dimensions().value();
            for (const auto& dimension : dimensions) {
                domain.add_dimension(dimension);
            }
        }

        if (column->tiledb_enumerations().has_value()) {
            auto enumerations = column->tiledb_enumerations().value();
            for (const auto& enumeration : enumerations) {
                ArraySchemaExperimental::add_enumeration(*ctx, schema, enumeration);
            }
        }

        if (column->tiledb_attributes().has_value()) {
            auto attributes = column->tiledb_attributes().value();
            for (const auto& attribute : attributes) {
                schema.add_attribute(attribute);
            }
        }
    }

    for (const auto& column : columns | std::views::filter([](const auto& col) { return !col->isIndexColumn(); })) {
        column->serialize(soma_schema_extension[TILEDB_SOMA_SCHEMA_COL_KEY]);

        if (column->tiledb_enumerations().has_value()) {
            auto enumerations = column->tiledb_enumerations().value();
            for (const auto& enumeration : enumerations) {
                ArraySchemaExperimental::add_enumeration(*ctx, schema, enumeration);
            }
        }

        if (column->tiledb_attributes().has_value()) {
            auto attributes = column->tiledb_attributes().value();
            for (const auto& attribute : attributes) {
                schema.add_attribute(attribute);
            }
        }
    }

    LOG_DEBUG(fmt::format("[ArrowAdapter] set_domain"));
    schema.set_domain(domain);

    LOG_DEBUG(fmt::format("[ArrowAdapter] index_column_info length {}", index_column_array->length));

    // Note: this must be done after we've got the core domain, since the
    // NDRectangle constructor requires access to the core domain.

    CurrentDomain current_domain(*ctx);
    NDRectangle ndrect(*ctx, domain);

    for (auto column : columns) {
        if (!column->isIndexColumn()) {
            continue;
        }

        column->set_current_domain_slot(ndrect, get_table_any_column_by_name<2>(index_column_info, column->name(), 3));

        // if (column->name() == SOMA_GEOMETRY_COLUMN_NAME) {
        //     std::vector<std::any> cdslot;
        //     for (int64_t j = 0; j < spatial_column_info.first->n_children;
        //          ++j) {
        //         cdslot.push_back(ArrowAdapter::get_table_any_column<2>(
        //             spatial_column_info.first->children[j],
        //             spatial_column_info.second->children[j],
        //             3));
        //     }

        //     column->set_current_domain_slot(ndrect, cdslot);
        // } else {
        //     column->set_current_domain_slot(
        //         ndrect,
        //         get_table_any_column_by_name<2>(
        //             index_column_info, column->name(), 3));
        // }
    }
    current_domain.set_ndrectangle(ndrect);

    LOG_DEBUG(fmt::format("[ArrowAdapter] before setting current_domain from ndrect"));
    ArraySchemaExperimental::set_current_domain(*ctx, schema, current_domain);
    LOG_DEBUG(fmt::format("[ArrowAdapter] after setting current_domain from ndrect"));

    LOG_DEBUG(fmt::format("[ArrowAdapter] check"));
    schema.check();

    LOG_DEBUG(fmt::format("[ArrowAdapter] returning"));
    return std::make_tuple(schema, soma_schema_extension);
}

std::pair<Attribute, std::optional<Enumeration>> ArrowAdapter::tiledb_attribute_from_arrow_schema(
    std::shared_ptr<Context> ctx,
    ArrowSchema* arrow_schema,
    std::string_view type_metadata,
    PlatformConfig platform_config) {
    auto type = ArrowAdapter::to_tiledb_format(arrow_schema->format, type_metadata);

    Attribute attr(*ctx, arrow_schema->name, type);

    FilterList filter_list = utils::create_attr_filter_list(arrow_schema->name, platform_config, ctx);
    attr.set_filter_list(filter_list);

    if (arrow_schema->flags & ARROW_FLAG_NULLABLE) {
        attr.set_nullable(true);
    }

    if (ArrowAdapter::arrow_is_var_length_type(arrow_schema->format)) {
        attr.set_cell_val_num(TILEDB_VAR_NUM);
    }

    std::optional<Enumeration> enmr = std::nullopt;

    if (arrow_schema->dictionary != nullptr) {
        auto enmr_format = arrow_schema->dictionary->format;
        auto enmr_type = ArrowAdapter::to_tiledb_format(enmr_format);
        auto enmr_label = util::get_enmr_label(arrow_schema, arrow_schema->dictionary);
        enmr = Enumeration::create_empty(
            *ctx,
            enmr_label,
            enmr_type,
            ArrowAdapter::arrow_is_var_length_type(enmr_format) ? TILEDB_VAR_NUM : 1,
            arrow_schema->flags & ARROW_FLAG_DICTIONARY_ORDERED);
        AttributeExperimental::set_enumeration_name(*ctx, attr, enmr_label);
        LOG_DEBUG(
            fmt::format(
                "[ArrowAdapter] dictionary for '{}' as '{}' '{}'",
                std::string(arrow_schema->name),
                tiledb::impl::type_to_str(enmr_type),
                std::string(enmr_format)));
    }

    return {attr, enmr};
}

inline void exitIfError(const ArrowErrorCode ec, const std::string& msg) {
    if (ec != NANOARROW_OK)
        throw TileDBSOMAError(fmt::format("ArrowAdapter: Arrow Error {} ", msg));
}

std::pair<managed_unique_ptr<ArrowArray>, managed_unique_ptr<ArrowSchema>> ArrowAdapter::to_arrow(
    std::shared_ptr<ReadColumnBuffer> column, bool downcast_dict_of_large_var) {
    managed_unique_ptr<ArrowSchema> schema = make_managed_unique<ArrowSchema>();
    managed_unique_ptr<ArrowArray> array = make_managed_unique<ArrowArray>();
    auto sch = schema.get();
    auto arr = array.get();

    auto coltype = to_arrow_format(column->type(), true).data();
    auto natype = to_nanoarrow_type(coltype);

    std::unordered_map<std::string, std::shared_ptr<ArrowBuffer>> enmr_map;

    if (natype == NANOARROW_TYPE_TIMESTAMP) {
        ArrowSchemaInit(sch);
        auto [ts_type, ts_unit] = to_nanoarrow_time(coltype);
        exitIfError(ArrowSchemaSetTypeDateTime(sch, ts_type, ts_unit, NULL), "Bad datetime");
    } else {
        exitIfError(ArrowSchemaInitFromType(sch, natype), "Bad schema init");
    }

    exitIfError(ArrowSchemaSetName(sch, column->name().data()), "Bad schema name");
    exitIfError(ArrowSchemaAllocateChildren(sch, 0), "Bad schema children alloc");
    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    schema->release = &release_schema;

    // this will be 3 for char vecs and 2 for enumerations
    int n_buffers = column->is_var() ? 3 : 2;

    exitIfError(ArrowArrayInitFromSchema(arr, sch, NULL), "Bad array init");
    exitIfError(ArrowArrayAllocateChildren(arr, 0), "Bad array children alloc");
    array->length = column->size();

    if (column->is_nullable()) {
        schema->flags |= ARROW_FLAG_NULLABLE;

        // Count nulls
        for (size_t i = 0; i < column->size(); ++i) {
            array->null_count += column->validity()[i] == 0;
        }
    } else {
        schema->flags &= ~ARROW_FLAG_NULLABLE;
    }

    // Create an ArrowBuffer to manage the lifetime of `column`.
    // - `arrow_buffer` holds shared_ptr to `column`, increments
    //   the use count and keeps the ColumnBuffer data alive.
    // - When the arrow array is released, `array->release()` is
    //   called with `arrow_buffer` in `private_data`.
    //   `arrow_buffer` is deleted, which decrements the the
    //   `column` use count. When the `column` use count reaches
    //   0, the ColumnBuffer data will be deleted.
    auto arrow_buffer = new PrivateArrowBuffer(std::make_shared<ArrowBuffer>(column.get()));

    LOG_TRACE(
        fmt::format(
            "[ArrowAdapter] column type {} name {} nbuf {} {} nullable {}",
            to_arrow_format(column->type()).data(),
            column->name().data(),
            n_buffers,
            array->n_buffers,
            column->is_nullable()));

    if (array->n_buffers != n_buffers) {
        throw TileDBSOMAError(
            fmt::format(
                "[ArrowAdapter] expected array n_buffers {} for column {}; got {}",
                n_buffers,
                column->name(),
                array->n_buffers));
    }

    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    array->release = &release_array;
    if (array->private_data != nullptr) {  // as we use nanoarrow's init
        free(array->private_data);         // free what was allocated before
    }  // assigning our ArrowBuffer pointer
    array->private_data = (void*)arrow_buffer;

    LOG_TRACE(fmt::format("[ArrowAdapter] create array name='{}' use_count={}", column->name(), column.use_count()));

    array->buffers = (const void**)malloc(sizeof(void*) * n_buffers);
    assert(array->buffers != nullptr);
    array->buffers[0] = nullptr;  // validity addressed below
    array->buffers[n_buffers - 1] = arrow_buffer->buffer_->data_.get();
    if (n_buffers == 3) {
        array->buffers[1] = arrow_buffer->buffer_->large_offsets_.get();
    }

    if (column->is_nullable()) {
        array->buffers[0] = arrow_buffer->buffer_->validity_.get();
    }

    if (column->is_ordered()) {
        schema->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
    }

    auto enmr = column->get_enumeration_info();
    if (enmr.has_value()) {
        if (!enmr_map.contains(enmr->name())) {
            enmr_map.insert(
                std::make_pair(enmr->name(), std::make_shared<ArrowBuffer>(*enmr, !downcast_dict_of_large_var)));
        }

        PrivateArrowBuffer* enmr_buffer = new PrivateArrowBuffer(enmr_map.at(enmr->name()));
        auto dict_sch = (ArrowSchema*)malloc(sizeof(ArrowSchema));
        auto dict_arr = (ArrowArray*)malloc(sizeof(ArrowArray));

        auto dcoltype = to_arrow_format(enmr->type(), !downcast_dict_of_large_var).data();
        auto dnatype = to_nanoarrow_type(dcoltype);

        if (dnatype == NANOARROW_TYPE_TIMESTAMP) {
            ArrowSchemaInit(dict_sch);
            auto [ts_type, ts_unit] = to_nanoarrow_time(dcoltype);
            exitIfError(ArrowSchemaSetTypeDateTime(dict_sch, ts_type, ts_unit, NULL), "Bad datetime");
        } else {
            exitIfError(ArrowSchemaInitFromType(dict_sch, dnatype), "Bad dict schema init");
        }

        exitIfError(ArrowSchemaSetName(dict_sch, ""), "Bad dict schema name");
        exitIfError(ArrowSchemaAllocateChildren(dict_sch, 0), "Bad dict schema children alloc");
        dict_sch->release = &release_schema;

        exitIfError(ArrowArrayInitFromSchema(dict_arr, dict_sch, NULL), "Bad dict array init");
        exitIfError(ArrowArrayAllocateChildren(dict_arr, 0), "Bad array children alloc");
        dict_arr->release = &release_array;

        bool is_var_enum = enmr->type() == TILEDB_STRING_ASCII || enmr->type() == TILEDB_STRING_UTF8 ||
                           enmr->type() == TILEDB_CHAR || enmr->type() == TILEDB_BLOB;
        int n_buffers_enum = is_var_enum ? 3 : 2;

        dict_arr->n_buffers = n_buffers_enum;
        dict_arr->buffers = (const void**)malloc(sizeof(void*) * n_buffers_enum);
        assert(dict_arr->buffers != nullptr);

        dict_arr->buffers[0] = nullptr;
        dict_arr->buffers[dict_arr->n_buffers - 1] = enmr_buffer->buffer_->data_.get();
        if (is_var_enum) {
            if (downcast_dict_of_large_var) {
                dict_arr->buffers[1] = enmr_buffer->buffer_->small_offsets_.get();
            } else {
                dict_arr->buffers[1] = enmr_buffer->buffer_->large_offsets_.get();
            }
        }

        dict_arr->length = enmr_buffer->buffer_->length;
        if (dict_arr->private_data != nullptr) {  // as we use nanoarrow's init
            free(dict_arr->private_data);         // free what was allocated before
        }  // assigning our ArrowBuffer pointer
        dict_arr->private_data = (void*)enmr_buffer;

        schema->dictionary = dict_sch;
        array->dictionary = dict_arr;
    }

    return std::pair(std::move(array), std::move(schema));
}

bool ArrowAdapter::arrow_is_var_length_type(const char* format) {
    return (
        (strcmp(format, "U") == 0) || (strcmp(format, "Z") == 0) || (strcmp(format, "u") == 0) ||
        (strcmp(format, "z") == 0));
}

std::string_view ArrowAdapter::to_arrow_format(tiledb_datatype_t tiledb_dtype, bool use_large) {
    auto u = use_large ? "U" : "u";
    auto z = use_large ? "Z" : "z";
    std::map<tiledb_datatype_t, std::string_view> _to_arrow_format_map = {
        {TILEDB_STRING_ASCII, u},     {TILEDB_CHAR, z},
        {TILEDB_STRING_UTF8, u},      {TILEDB_BLOB, z},
        {TILEDB_INT8, "c"},           {TILEDB_UINT8, "C"},
        {TILEDB_INT16, "s"},          {TILEDB_UINT16, "S"},
        {TILEDB_INT32, "i"},          {TILEDB_UINT32, "I"},
        {TILEDB_INT64, "l"},          {TILEDB_UINT64, "L"},
        {TILEDB_FLOAT32, "f"},        {TILEDB_FLOAT64, "g"},
        {TILEDB_BOOL, "b"},           {TILEDB_DATETIME_SEC, "tss:"},
        {TILEDB_DATETIME_MS, "tsm:"}, {TILEDB_DATETIME_US, "tsu:"},
        {TILEDB_DATETIME_NS, "tsn:"}, {TILEDB_GEOM_WKB, z},
        {TILEDB_GEOM_WKT, u}};
    try {
        return _to_arrow_format_map.at(tiledb_dtype);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(
            fmt::format("ArrowAdapter: Unsupported TileDB type: {} ", tiledb::impl::type_to_str(tiledb_dtype)));
    }
}

tiledb_datatype_t ArrowAdapter::to_tiledb_format(std::string_view arrow_dtype, std::string_view arrow_dtype_metadata) {
    std::map<std::string_view, tiledb_datatype_t> _to_tiledb_format_map = {
        {"u", TILEDB_STRING_UTF8},    {"U", TILEDB_STRING_UTF8},
        {"z", TILEDB_CHAR},           {"Z", TILEDB_CHAR},
        {"c", TILEDB_INT8},           {"C", TILEDB_UINT8},
        {"s", TILEDB_INT16},          {"S", TILEDB_UINT16},
        {"i", TILEDB_INT32},          {"I", TILEDB_UINT32},
        {"l", TILEDB_INT64},          {"L", TILEDB_UINT64},
        {"f", TILEDB_FLOAT32},        {"g", TILEDB_FLOAT64},
        {"b", TILEDB_BOOL},           {"tss:", TILEDB_DATETIME_SEC},
        {"tsm:", TILEDB_DATETIME_MS}, {"tsu:", TILEDB_DATETIME_US},
        {"tsn:", TILEDB_DATETIME_NS},
    };

    try {
        auto dtype = _to_tiledb_format_map.at(arrow_dtype);

        if (dtype == TILEDB_CHAR && arrow_dtype_metadata.compare("WKB") == 0) {
            dtype = TILEDB_GEOM_WKB;
        } else if (dtype == TILEDB_STRING_UTF8 && arrow_dtype_metadata.compare("WKT") == 0) {
            dtype = TILEDB_GEOM_WKT;
        }

        return dtype;
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(
            fmt::format("ArrowAdapter: Unsupported Arrow type: {} ({})", arrow_dtype, to_arrow_readable(arrow_dtype)));
    }
}

managed_unique_ptr<ArrowSchema> ArrowAdapter::make_arrow_schema(
    const std::vector<std::string>& names, const std::vector<tiledb_datatype_t>& tiledb_datatypes) {
    auto num_names = names.size();
    auto num_types = tiledb_datatypes.size();

    if (num_names != num_types) {
        throw TileDBSOMAError(
            fmt::format(
                "ArrowAdapter::make_arrow_schema: internal coding error: num_types "
                "{} != num_names {}",
                num_names,
                num_types));
    }

    auto arrow_schema = make_managed_unique<ArrowSchema>();
    arrow_schema->format = strdup("+s");  // structure, i.e. non-leaf node
    arrow_schema->name = strdup("parent");
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    arrow_schema->n_children = num_names;  // non-leaf node
    arrow_schema->children = (ArrowSchema**)malloc(arrow_schema->n_children * sizeof(ArrowSchema*));
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;

    LOG_DEBUG(fmt::format("[ArrowAdapter] make_arrow_schema n_children {}", arrow_schema->n_children));

    for (int i = 0; i < (int)num_names; i++) {
        auto dim_schema = make_arrow_schema_child(names[i], tiledb_datatypes[i]);
        LOG_TRACE(
            fmt::format(
                "[ArrowAdapter] make_arrow_schema child {} format {} name {}",
                i,
                dim_schema->format,
                dim_schema->name));
        arrow_schema->children[i] = dim_schema;
    }

    return arrow_schema;
}

ArrowSchema* ArrowAdapter::make_arrow_schema_child(std::string name, tiledb_datatype_t tiledb_datatype) {
    ArrowSchema* arrow_schema = (ArrowSchema*)malloc(sizeof(ArrowSchema));
    auto arrow_type_name = ArrowAdapter::tdb_to_arrow_type(tiledb_datatype);
    arrow_schema->name = strdup(name.c_str());
    arrow_schema->format = strdup(arrow_type_name.c_str());
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    arrow_schema->n_children = 0;      // leaf node
    arrow_schema->children = nullptr;  // leaf node
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;

    if (strcmp(arrow_schema->name, SOMA_GEOMETRY_COLUMN_NAME.c_str()) == 0) {
        nanoarrow::UniqueBuffer buffer;
        ArrowMetadataBuilderInit(buffer.get(), nullptr);
        ArrowMetadataBuilderAppend(
            buffer.get(),
            ArrowCharView(ARROW_DATATYPE_METADATA_KEY.c_str()),
            ArrowCharView(tiledb_datatype == TILEDB_GEOM_WKB ? "WKB" : "WKT"));
        ArrowSchemaSetMetadata(arrow_schema, std::string((char*)buffer->data, buffer->size_bytes).c_str());
    }

    return arrow_schema;
}

managed_unique_ptr<ArrowSchema> ArrowAdapter::make_arrow_schema_parent(size_t num_columns, std::string_view name) {
    auto arrow_schema = make_managed_unique<ArrowSchema>();

    arrow_schema->format = strdup("+s");  // structure, i.e. non-leaf node
    arrow_schema->name = strdup(name.data());
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    arrow_schema->n_children = static_cast<int64_t>(num_columns);  // non-leaf node
    arrow_schema->children = (ArrowSchema**)malloc(arrow_schema->n_children * sizeof(ArrowSchema*));
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;

    for (size_t i = 0; i < num_columns; i++) {
        arrow_schema->children[i] = nullptr;
    }

    LOG_DEBUG(fmt::format("[ArrowAdapter] make_arrow_schema n_children {}", arrow_schema->n_children));

    return arrow_schema;
}

managed_unique_ptr<ArrowArray> ArrowAdapter::make_arrow_array_parent(size_t num_columns) {
    auto arrow_array = make_managed_unique<ArrowArray>();

    // All zero/null since this is a parent ArrowArray, and each
    // column/child is also of type ArrowArray.
    arrow_array->length = 0;
    arrow_array->null_count = 0;
    arrow_array->offset = 0;
    arrow_array->n_buffers = 0;
    arrow_array->n_children = static_cast<int64_t>(num_columns);
    arrow_array->buffers = nullptr;
    arrow_array->dictionary = nullptr;
    arrow_array->release = &ArrowAdapter::release_array;
    arrow_array->private_data = nullptr;

    arrow_array->children = (ArrowArray**)malloc(num_columns * sizeof(ArrowArray*));
    for (size_t i = 0; i < num_columns; i++) {
        arrow_array->children[i] = nullptr;
    }

    LOG_DEBUG(fmt::format("[ArrowAdapter] make_arrow_array n_children {}", arrow_array->n_children));

    return arrow_array;
}

void ArrowAdapter::log_make_arrow_array_child(ArrowArray* child) {
    LOG_TRACE(
        fmt::format("[ArrowAdapter] make_arrow_array_child length {} n_buffers {}", child->length, child->n_buffers));
}

void ArrowAdapter::_check_shapes(ArrowArray* arrow_array, ArrowSchema* arrow_schema) {
    if (arrow_array->n_children != arrow_schema->n_children) {
        throw std::runtime_error(
            "ArrowAdapter::_check_shapes: internal coding error: data/schema "
            "mismatch");
    }
    for (int64_t i = 0; i < arrow_array->n_children; i++) {
        _check_shapes(arrow_array->children[i], arrow_schema->children[i]);
    }
}

int64_t ArrowAdapter::_get_column_index_from_name(const ArrowTable& arrow_table, std::string column_name) {
    ArrowArray* arrow_array = arrow_table.first.get();
    ArrowSchema* arrow_schema = arrow_table.second.get();
    // Make sure the child-count is the same
    _check_shapes(arrow_array, arrow_schema);

    if (arrow_schema->n_children == 0) {
        throw std::runtime_error(
            "ArrowAdapter::_check_shapes: internal coding error: childless "
            "table");
    }

    for (int64_t i = 0; i < arrow_schema->n_children; i++) {
        if (strcmp(arrow_schema->children[i]->name, column_name.c_str()) == 0) {
            return i;
        }
    }

    throw std::runtime_error(fmt::format("ArrowAdapter::_check_shapes: column {} not found", column_name));
}

ArrowArray* ArrowAdapter::_get_and_check_column(
    const ArrowTable& arrow_table, int64_t column_index, int64_t expected_n_buffers) {
    ArrowArray* arrow_array = arrow_table.first.get();
    if (column_index < 0 || column_index >= arrow_array->n_children) {
        throw std::runtime_error(
            fmt::format(
                "ArrowAdapter::_get_and_check_column: column index {} out of "
                "bounds {}..{}",
                column_index,
                0,
                arrow_array->n_children - 1));
    }

    ArrowArray* child = arrow_array->children[column_index];

    if (child->n_children != 0) {
        throw std::runtime_error(
            fmt::format(
                "ArrowAdapter::_get_and_check_column: column index {} is "
                "non-terminal",
                column_index));
    }

    if (expected_n_buffers == 2) {
        if (child->n_buffers != 2) {
            throw std::runtime_error(
                fmt::format(
                    "ArrowAdapter::_get_and_check_column: column index {} "
                    "has buffer count {}; expected 2 for non-string data",
                    column_index,
                    child->n_buffers));
        }

    } else if (expected_n_buffers == 3) {
        if (child->n_buffers != 3) {
            throw std::runtime_error(
                fmt::format(
                    "ArrowAdapter::_get_and_check_column: column index {} is "
                    "has buffer count {}; expected 3 for string data",
                    column_index,
                    child->n_buffers));
        }

    } else {
        throw std::runtime_error(
            fmt::format(
                "ArrowAdapter::_get_and_check_column: internal coding error: "
                "expected_n_buffers {} is "
                "neither 2 nor 3.",
                expected_n_buffers));
    }

    return child;
}

managed_unique_ptr<ArrowArray> ArrowAdapter::arrow_array_insert_at_index(
    managed_unique_ptr<ArrowArray> parent_array,
    std::vector<managed_unique_ptr<ArrowArray>> child_arrays,
    int64_t index) {
    if (parent_array->n_children < index || index < 0) {
        throw std::runtime_error(
            "[ArrowAdapter][arrow_array_insert_at_index] Invalid index to "
            "insert array");
    }

    if (child_arrays.size() == 0) {
        return parent_array;
    }

    auto array = make_arrow_array_parent(static_cast<size_t>(parent_array->n_children) + child_arrays.size());

    for (int64_t i = 0; i < array->n_children; ++i) {
        int64_t idx = i <= index ? i : i - child_arrays.size();
        array->children[i] = (ArrowArray*)malloc(sizeof(ArrowArray));

        if (i >= index && i < index + static_cast<int64_t>(child_arrays.size())) {
            ArrowArrayMove(child_arrays[i - index].get(), array->children[i]);
        } else {
            ArrowArrayMove(parent_array->children[idx], array->children[i]);
        }
    }

    return array;
}

managed_unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_insert_at_index(
    managed_unique_ptr<ArrowSchema> parent_schema,
    std::vector<managed_unique_ptr<ArrowSchema>> child_schemas,
    int64_t index) {
    if (parent_schema->n_children < index || index < 0) {
        throw std::runtime_error(
            "[ArrowAdapter][arrow_schema_insert_at_index] Invalid index to "
            "insert schema");
    }

    if (child_schemas.size() == 0) {
        return parent_schema;
    }

    auto schema = make_arrow_schema_parent(parent_schema->n_children + child_schemas.size());

    for (int64_t i = 0; i < schema->n_children; ++i) {
        int64_t idx = i <= index ? i : i - child_schemas.size();
        schema->children[i] = (ArrowSchema*)malloc(sizeof(ArrowSchema));

        if (i >= index && i < index + static_cast<int64_t>(child_schemas.size())) {
            ArrowSchemaMove(child_schemas[i - index].get(), schema->children[i]);
        } else {
            ArrowSchemaMove(parent_schema->children[idx], schema->children[i]);
        }
    }

    return schema;
}

managed_unique_ptr<ArrowArray> ArrowAdapter::arrow_array_remove_at_index(
    managed_unique_ptr<ArrowArray> array, int64_t index) {
    if (array->n_children <= index || index < 0) {
        throw std::runtime_error(
            "[ArrowAdapter][arrow_array_remove_at_index] Invalid index to "
            "remove child array");
    }

    auto array_new = make_arrow_array_parent(array->n_children - 1);
    for (int64_t i = 0; i < array->n_children; ++i) {
        int64_t idx = i <= index ? i : i - 1;

        if (i != index) {
            array_new->children[idx] = (ArrowArray*)malloc(sizeof(ArrowArray));
            ArrowArrayMove(array->children[i], array_new->children[idx]);
        }
    }

    return array_new;
}

managed_unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_remove_at_index(
    managed_unique_ptr<ArrowSchema> schema, int64_t index) {
    if (schema->n_children <= index || index < 0) {
        throw std::runtime_error(
            "[ArrowAdapter][arrow_schema_remove_at_index] Invalid index to "
            "remove child schema");
    }

    auto schema_new = make_arrow_schema_parent(schema->n_children - 1);

    for (int64_t i = 0; i < schema->n_children; ++i) {
        int64_t idx = i <= index ? i : i - 1;

        if (i != index) {
            schema_new->children[idx] = (ArrowSchema*)malloc(sizeof(ArrowSchema));
            ArrowSchemaMove(schema->children[i], schema_new->children[idx]);
        }
    }

    return schema_new;
}

size_t ArrowAdapter::_set_var_dictionary_buffers(
    Enumeration& enumeration, const Context& ctx, const void** buffers, bool downcast_dict_of_large_var) {
    const void* data;
    uint64_t data_size;

    ctx.handle_error(tiledb_enumeration_get_data(ctx.ptr().get(), enumeration.ptr().get(), &data, &data_size));

    const void* offsets;
    uint64_t offsets_size;
    ctx.handle_error(tiledb_enumeration_get_offsets(ctx.ptr().get(), enumeration.ptr().get(), &offsets, &offsets_size));

    size_t count = offsets_size / sizeof(uint64_t);

    if (downcast_dict_of_large_var) {
        std::span<const uint64_t> offsets_v(static_cast<const uint64_t*>(offsets), count);

        uint32_t* small_offsets = static_cast<uint32_t*>(malloc((count + 1) * sizeof(uint32_t)));
        for (size_t i = 0; i < count; ++i) {
            small_offsets[i] = static_cast<uint32_t>(offsets_v[i]);
        }
        small_offsets[count] = static_cast<uint32_t>(data_size);
        buffers[1] = small_offsets;
    } else {
        auto large_offsets = static_cast<uint64_t*>(malloc((count + 1) * sizeof(uint64_t)));
        std::memcpy(static_cast<void*>(large_offsets), offsets, offsets_size);
        large_offsets[count] = data_size;
        buffers[1] = large_offsets;
    }

    buffers[2] = malloc(data_size);
    std::memcpy(const_cast<void*>(buffers[2]), data, data_size);
    return count;
}

size_t ArrowAdapter::_set_dictionary_buffers(Enumeration& enumeration, const Context& ctx, const void** buffers) {
    const void* data;
    uint64_t data_size;

    ctx.handle_error(tiledb_enumeration_get_data(ctx.ptr().get(), enumeration.ptr().get(), &data, &data_size));

    buffers[1] = malloc(data_size);
    std::memcpy(const_cast<void*>(buffers[1]), data, data_size);

    switch (enumeration.type()) {
        case TILEDB_INT8:
            return data_size / sizeof(int8_t);
        case TILEDB_UINT8:
            return data_size / sizeof(uint8_t);
        case TILEDB_INT16:
            return data_size / sizeof(int16_t);
        case TILEDB_UINT16:
            return data_size / sizeof(uint16_t);
        case TILEDB_INT32:
            return data_size / sizeof(int32_t);
        case TILEDB_UINT32:
            return data_size / sizeof(uint32_t);
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
        case TILEDB_INT64:
            return data_size / sizeof(int64_t);
        case TILEDB_UINT64:
            return data_size / sizeof(uint64_t);
        case TILEDB_FLOAT32:
            return data_size / sizeof(float_t);
        case TILEDB_FLOAT64:
            return data_size / sizeof(double_t);
        default:
            throw TileDBSOMAError(
                fmt::format(
                    "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                    tiledb::impl::type_to_str(enumeration.type())));
    }
}

size_t ArrowAdapter::_set_bool_dictionary_buffers(Enumeration& enumeration, const Context& ctx, const void** buffers) {
    const void* data;
    uint64_t data_size;

    ctx.handle_error(tiledb_enumeration_get_data(ctx.ptr().get(), enumeration.ptr().get(), &data, &data_size));

    std::span<const bool> data_v(static_cast<const bool*>(data), data_size);
    size_t count = data_size / sizeof(bool);

    // Represent the Boolean vector with, at most, the last two
    // bits. In Arrow, Boolean values are LSB packed
    uint8_t packed_data = 0;
    for (size_t i = 0; i < count; ++i)
        packed_data |= (data_v[i] << i);

    // Allocate a single byte to copy the bits into
    buffers[1] = malloc(1);
    std::memcpy(const_cast<void*>(buffers[1]), &packed_data, 1);

    return count;
}

}  // namespace tiledbsoma
