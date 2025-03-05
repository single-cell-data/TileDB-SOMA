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

void ArrowAdapter::release_schema(struct ArrowSchema* schema) {
    std::string name_for_log(
        schema->name == nullptr ? "anonymous" : schema->name);
    if (schema->name != nullptr)
        LOG_DEBUG(std::format(
            "[ArrowAdapter] release_schema start for {}", schema->name));

    if (schema->name != nullptr) {
        LOG_TRACE(std::format(
            "[ArrowAdapter] release_schema schema->name {}", schema->name));
        free((void*)schema->name);
        schema->name = nullptr;
    }
    if (schema->format != nullptr) {
        LOG_TRACE(std::format(
            "[ArrowAdapter] release_schema name {} schema->format {}",
            name_for_log,
            schema->format));
        free((void*)schema->format);
        schema->format = nullptr;
    }
    if (schema->metadata != nullptr) {
        LOG_TRACE(std::format(
            "[ArrowAdapter] release_schema name {} schema->metadata",
            name_for_log));
        free((void*)schema->metadata);
        schema->metadata = nullptr;
    }

    if (schema->children != nullptr) {
        LOG_TRACE(std::format(
            "[ArrowAdapter] release_schema name {} n_children {} begin "
            "recurse ",
            name_for_log,
            schema->n_children));

        for (auto i = 0; i < schema->n_children; i++) {
            if (schema->children[i] != nullptr) {
                if (schema->children[i]->release != nullptr) {
                    LOG_TRACE(std::format(
                        "[ArrowAdapter] release_schema name {} schema->child "
                        "{} "
                        "release",
                        name_for_log,
                        i));
                    schema->children[i]->release(schema->children[i]);
                }
                LOG_TRACE(std::format(
                    "[ArrowAdapter] release_schema name {} schema->child {} "
                    "free",
                    name_for_log,
                    i));
                free(schema->children[i]);
                schema->children[i] = nullptr;
            }
        }

        LOG_TRACE(std::format(
            "[ArrowAdapter] release_schema name {} n_children {} end recurse ",
            name_for_log,
            schema->n_children));

        free(schema->children);
        schema->children = nullptr;
    }

    if (schema->dictionary != nullptr) {
        if (schema->dictionary->release != nullptr) {
            LOG_TRACE(std::format(
                "[ArrowAdapter] release_schema name {} schema->dict release",
                name_for_log));
            release_schema(schema->dictionary);
        }
        LOG_TRACE(std::format(
            "[ArrowAdapter] release_schema name {} schema->dict free",
            name_for_log));
        free(schema->dictionary);
        schema->dictionary = nullptr;
    }

    schema->release = nullptr;
    LOG_TRACE(std::format(
        "[ArrowAdapter] release_schema name {} done", name_for_log));
}

void ArrowAdapter::release_array(struct ArrowArray* array) {
    auto arrow_buffer = static_cast<ArrowBuffer*>(array->private_data);
    if (arrow_buffer != nullptr) {
        LOG_TRACE(std::format(
            "[ArrowAdapter] release_array {} use_count={}",
            arrow_buffer->buffer_->name(),
            arrow_buffer->buffer_.use_count()));

        // Delete the ArrowBuffer, which was allocated with new.
        // If the ArrowBuffer.buffer_ shared_ptr is the last reference to the
        // underlying ColumnBuffer, the ColumnBuffer will be deleted.
        delete arrow_buffer;
    }

    if (array->buffers != nullptr) {
        free(array->buffers);
        array->buffers = nullptr;
    }

    if (array->children != nullptr) {
        for (auto i = 0; i < array->n_children; i++) {
            if (array->children[i] != nullptr) {
                if (array->children[i]->release != nullptr) {
                    LOG_TRACE(std::format(
                        "[ArrowAdapter] release_schema array->child {} release",
                        i));

                    array->children[i]->release(array->children[i]);
                }
                LOG_TRACE(std::format(
                    "[ArrowAdapter] release_schema array->child {} free", i));
                free(array->children[i]);
                array->children[i] = nullptr;
            }
        }
        LOG_TRACE("[ArrowAdapter] release_array array->children");
        free(array->children);
        array->children = nullptr;
    }

    if (array->dictionary != nullptr) {
        // Dictionary arrays are allocated differently than data arrays
        // We need to free the buffers one at a time, then we can call the
        // release schema to continue the cleanup properly
        for (size_t i = 0; i < array->dictionary->n_buffers; ++i) {
            if (array->dictionary->buffers[i] != nullptr) {
                free(const_cast<void*>(array->dictionary->buffers[i]));
                array->dictionary->buffers[i] = nullptr;
            }
        }

        LOG_TRACE("[ArrowAdapter] release_array array->dict release");

        array->dictionary->release(array->dictionary);
        free(array->dictionary);
        array->dictionary = nullptr;
    }

    array->release = nullptr;
    LOG_TRACE(std::format("[ArrowAdapter] release_array done"));
}

PlatformConfig ArrowAdapter::platform_config_from_tiledb_schema(
    ArraySchema tiledb_schema) {
    std::map<tiledb_layout_t, std::string> layout_as_string{
        {TILEDB_ROW_MAJOR, "row-major"},
        {TILEDB_COL_MAJOR, "column-major"},
        {TILEDB_HILBERT, "hilbert"},
        {TILEDB_UNORDERED, "unordered"},
    };

    PlatformConfig platform_config;
    platform_config.capacity = tiledb_schema.capacity();
    platform_config.allows_duplicates = tiledb_schema.allows_dups();
    platform_config.tile_order = layout_as_string[tiledb_schema.tile_order()];
    platform_config.cell_order = layout_as_string[tiledb_schema.cell_order()];
    platform_config.offsets_filters = ArrowAdapter::_get_filter_list_json(
                                          tiledb_schema.offsets_filter_list())
                                          .dump();
    platform_config.validity_filters = ArrowAdapter::_get_filter_list_json(
                                           tiledb_schema.validity_filter_list())
                                           .dump();
    platform_config.attrs = ArrowAdapter::_get_attrs_filter_list_json(
                                tiledb_schema)
                                .dump();
    platform_config.dims = ArrowAdapter::_get_dims_list_json(tiledb_schema)
                               .dump();

    return platform_config;
}

PlatformSchemaConfig ArrowAdapter::platform_schema_config_from_tiledb(
    ArraySchema tiledb_schema) {
    std::map<tiledb_layout_t, std::string> layout_as_string{
        {TILEDB_ROW_MAJOR, "row-major"},
        {TILEDB_COL_MAJOR, "column-major"},
        {TILEDB_HILBERT, "hilbert"},
        {TILEDB_UNORDERED, "unordered"},
    };

    PlatformSchemaConfig platform_config;
    platform_config.capacity = tiledb_schema.capacity();
    platform_config.allows_duplicates = tiledb_schema.allows_dups();
    platform_config.tile_order = layout_as_string[tiledb_schema.tile_order()];
    platform_config.cell_order = layout_as_string[tiledb_schema.cell_order()];
    platform_config.offsets_filters = ArrowAdapter::_get_filter_list_json(
                                          tiledb_schema.offsets_filter_list())
                                          .dump();
    platform_config.validity_filters = ArrowAdapter::_get_filter_list_json(
                                           tiledb_schema.validity_filter_list())
                                           .dump();
    platform_config.attrs = ArrowAdapter::_get_attrs_filter_list_json(
                                tiledb_schema)
                                .dump();
    platform_config.dims = ArrowAdapter::_get_dims_list_json(tiledb_schema)
                               .dump();

    return platform_config;
}

json ArrowAdapter::_get_attrs_filter_list_json(
    const ArraySchema& tiledb_schema) {
    json attrs_filter_list_as_json;
    for (const auto& attr : tiledb_schema.attributes()) {
        json attr_info = {
            {"filters", _get_filter_list_json(attr.second.filter_list())}};
        attrs_filter_list_as_json.emplace(attr.first, attr_info);
    }
    return attrs_filter_list_as_json;
}

json ArrowAdapter::_get_dims_list_json(const ArraySchema& tiledb_schema) {
    json dims_as_json;
    for (const auto& dim : tiledb_schema.domain().dimensions()) {
        json dim_info = {
            {"tile", dim.tile_extent_to_str()},
            {"filters", _get_filter_list_json(dim.filter_list())}};
        dims_as_json.emplace(dim.name(), dim_info);
    }
    return dims_as_json;
}

json ArrowAdapter::_get_filter_list_json(FilterList filter_list) {
    std::map<tiledb_filter_option_t, std::string> option_as_string = {
        {TILEDB_COMPRESSION_LEVEL, "COMPRESSION_LEVEL"},
        {TILEDB_BIT_WIDTH_MAX_WINDOW, "BIT_WIDTH_MAX_WINDOW"},
        {TILEDB_POSITIVE_DELTA_MAX_WINDOW, "POSITIVE_DELTA_MAX_WINDOW"},
        {TILEDB_SCALE_FLOAT_BYTEWIDTH, "SCALE_FLOAT_BYTEWIDTH"},
        {TILEDB_SCALE_FLOAT_FACTOR, "SCALE_FLOAT_FACTOR"},
        {TILEDB_SCALE_FLOAT_OFFSET, "SCALE_FLOAT_OFFSET"},
        {TILEDB_WEBP_INPUT_FORMAT, "WEBP_INPUT_FORMAT"},
        {TILEDB_WEBP_QUALITY, "WEBP_QUALITY"},
        {TILEDB_WEBP_LOSSLESS, "WEBP_LOSSLESS"},
        {TILEDB_COMPRESSION_REINTERPRET_DATATYPE,
         "COMPRESSION_REINTERPRET_DATATYPE"},
    };

    json filter_list_as_json = {};
    for (uint32_t i = 0; i < filter_list.nfilters(); ++i) {
        json filter_as_json = {};

        auto filter = filter_list.filter(i);
        filter_as_json.emplace("name", Filter::to_str(filter.filter_type()));

        switch (filter.filter_type()) {
            case TILEDB_FILTER_GZIP:
            case TILEDB_FILTER_ZSTD:
            case TILEDB_FILTER_LZ4:
            case TILEDB_FILTER_BZIP2:
            case TILEDB_FILTER_RLE:
            case TILEDB_FILTER_DICTIONARY:
                filter_as_json.emplace(
                    "COMPRESSION_LEVEL",
                    filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL));
                break;

            case TILEDB_FILTER_DELTA:
            case TILEDB_FILTER_DOUBLE_DELTA:
                filter_as_json.emplace(
                    "COMPRESSION_LEVEL",
                    filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL));
                filter_as_json.emplace(
                    "COMPRESSION_REINTERPRET_DATATYPE",
                    filter.get_option<uint8_t>(
                        TILEDB_COMPRESSION_REINTERPRET_DATATYPE));
                break;

            case TILEDB_FILTER_BIT_WIDTH_REDUCTION:
                filter_as_json.emplace(
                    "BIT_WIDTH_MAX_WINDOW",
                    filter.get_option<uint32_t>(TILEDB_BIT_WIDTH_MAX_WINDOW));
                break;

            case TILEDB_FILTER_POSITIVE_DELTA:
                filter_as_json.emplace(
                    "POSITIVE_DELTA_MAX_WINDOW",
                    filter.get_option<uint32_t>(
                        TILEDB_POSITIVE_DELTA_MAX_WINDOW));
                break;

            case TILEDB_FILTER_SCALE_FLOAT:
                filter_as_json.emplace(
                    "SCALE_FLOAT_FACTOR",
                    filter.get_option<double>(TILEDB_SCALE_FLOAT_FACTOR));
                filter_as_json.emplace(
                    "SCALE_FLOAT_OFFSET",
                    filter.get_option<double>(TILEDB_SCALE_FLOAT_OFFSET));
                filter_as_json.emplace(
                    "SCALE_FLOAT_BYTEWIDTH",
                    filter.get_option<uint64_t>(TILEDB_SCALE_FLOAT_BYTEWIDTH));
                break;

            case TILEDB_FILTER_WEBP:
                filter_as_json.emplace(
                    "WEBP_INPUT_FORMAT",
                    filter.get_option<uint8_t>(TILEDB_WEBP_INPUT_FORMAT));
                filter_as_json.emplace(
                    "WEBP_QUALITY",
                    filter.get_option<float>(TILEDB_WEBP_QUALITY));
                filter_as_json.emplace(
                    "WEBP_LOSSLESS",
                    filter.get_option<uint8_t>(TILEDB_WEBP_LOSSLESS));
                break;

            case TILEDB_FILTER_CHECKSUM_MD5:
            case TILEDB_FILTER_CHECKSUM_SHA256:
            case TILEDB_FILTER_XOR:
            case TILEDB_FILTER_BITSHUFFLE:
            case TILEDB_FILTER_BYTESHUFFLE:
            case TILEDB_FILTER_DEPRECATED:
            case TILEDB_FILTER_NONE:
                // These filters have no options and are left empty
                // intentionally
                break;
        }
        filter_list_as_json.emplace_back(filter_as_json);
    }
    return filter_list_as_json;
}

std::unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_from_tiledb_array(
    std::shared_ptr<Context> ctx, std::shared_ptr<Array> tiledb_array) {
    auto tiledb_schema = tiledb_array->schema();
    auto ndim = tiledb_schema.domain().ndim();
    auto nattr = tiledb_schema.attribute_num();

    std::unique_ptr<ArrowSchema> arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = strdup("+s");
    arrow_schema->name = strdup("parent");
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    arrow_schema->n_children = ndim + nattr;
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;

    arrow_schema->children = (ArrowSchema**)malloc(
        arrow_schema->n_children * sizeof(ArrowSchema*));
    LOG_DEBUG(std::format(
        "[ArrowAdapter] arrow_schema_from_tiledb_array n_children {}",
        arrow_schema->n_children));

    ArrowSchema* child = nullptr;

    for (uint32_t i = 0; i < ndim; ++i) {
        auto dim = tiledb_schema.domain().dimension(i);
        child = arrow_schema->children[i] = (ArrowSchema*)malloc(
            sizeof(ArrowSchema));
        child->format = strdup(
            ArrowAdapter::to_arrow_format(dim.type()).data());
        child->name = strdup(dim.name().c_str());
        child->metadata = nullptr;
        child->flags = 0;
        child->n_children = 0;
        child->children = nullptr;
        child->dictionary = nullptr;
        child->release = &ArrowAdapter::release_schema;
        child->private_data = nullptr;
        LOG_TRACE(std::format(
            "[ArrowAdapter] arrow_schema_from_tiledb_array dim {} format {} "
            "name {}",
            i,
            child->format,
            child->name));
    }

    for (uint32_t i = 0; i < nattr; ++i) {
        auto attr = tiledb_schema.attribute(i);
        child = arrow_schema->children[ndim + i] = (ArrowSchema*)malloc(
            sizeof(ArrowSchema));
        child->format = strdup(
            ArrowAdapter::to_arrow_format(attr.type()).data());
        child->name = strdup(attr.name().c_str());
        child->metadata = nullptr;
        child->flags = 0;
        if (attr.nullable()) {
            child->flags |= ARROW_FLAG_NULLABLE;
        } else {
            child->flags &= ~ARROW_FLAG_NULLABLE;
        }
        child->n_children = 0;
        child->children = nullptr;
        child->dictionary = nullptr;
        child->release = &ArrowAdapter::release_schema;
        child->private_data = nullptr;

        LOG_TRACE(std::format(
            "[ArrowAdapter] arrow_schema_from_tiledb_array attr {} format {} "
            "name {}",
            i,
            child->format,
            child->name));

        auto enmr_name = AttributeExperimental::get_enumeration_name(
            *ctx, attr);
        if (enmr_name.has_value()) {
            auto enmr = ArrayExperimental::get_enumeration(
                *ctx, *tiledb_array, enmr_name.value());
            auto dict = (ArrowSchema*)malloc(sizeof(ArrowSchema));
            dict->format = strdup(
                ArrowAdapter::to_arrow_format(enmr.type(), false).data());
            if (enmr.type() == TILEDB_STRING_ASCII ||
                enmr.type() == TILEDB_CHAR) {
                dict->format = strdup("z");
            } else {
                dict->format = strdup(
                    ArrowAdapter::to_arrow_format(enmr.type(), false).data());
            }
            dict->name = strdup(enmr.name().c_str());
            dict->metadata = nullptr;
            if (enmr.ordered()) {
                child->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
            } else {
                child->flags &= ~ARROW_FLAG_DICTIONARY_ORDERED;
            }
            dict->n_children = 0;
            dict->children = nullptr;
            dict->dictionary = nullptr;
            dict->release = &ArrowAdapter::release_schema;
            dict->private_data = nullptr;
            child->dictionary = dict;
        }
        child->release = &ArrowAdapter::release_schema;
    }

    return arrow_schema;
}

std::unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_from_tiledb_dimension(
    const Dimension& dimension) {
    std::unique_ptr<ArrowSchema> arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = strdup(
        ArrowAdapter::to_arrow_format(dimension.type()).data());
    arrow_schema->name = strdup(dimension.name().c_str());
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    arrow_schema->n_children = 0;
    arrow_schema->children = nullptr;
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;
    LOG_TRACE(std::format(
        "[ArrowAdapter] arrow_schema_from_tiledb_dimension format {} "
        "name {}",
        arrow_schema->format,
        arrow_schema->name));

    return arrow_schema;
}

std::unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_from_tiledb_attribute(
    const Attribute& attribute, const Context& ctx, const Array& tiledb_array) {
    std::unique_ptr<ArrowSchema> arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = strdup(
        ArrowAdapter::to_arrow_format(attribute.type()).data());
    arrow_schema->name = strdup(attribute.name().c_str());
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
        ArrowMetadataBuilderAppend(
            metadata_buffer.get(),
            ArrowCharView("dtype"),
            ArrowCharView("WKB"));
        ArrowSchemaSetMetadata(
            arrow_schema.get(), reinterpret_cast<char*>(metadata_buffer->data));
    }

    LOG_TRACE(std::format(
        "[ArrowAdapter] arrow_schema_from_tiledb_array format {} "
        "name {}",
        arrow_schema->format,
        arrow_schema->name));
    // We shouldn;t have to cast constness away. Maybe missing const qualifier
    // from AttributeExperimental::get_enumeration_name
    auto enmr_name = AttributeExperimental::get_enumeration_name(
        ctx, const_cast<Attribute&>(attribute));
    if (enmr_name.has_value()) {
        auto enmr = ArrayExperimental::get_enumeration(
            ctx, tiledb_array, *enmr_name);
        auto dict = (ArrowSchema*)malloc(sizeof(ArrowSchema));
        dict->format = strdup(
            ArrowAdapter::to_arrow_format(enmr.type(), false).data());
        if (enmr.type() == TILEDB_STRING_ASCII || enmr.type() == TILEDB_CHAR) {
            dict->format = strdup("z");
        } else {
            dict->format = strdup(
                ArrowAdapter::to_arrow_format(enmr.type(), false).data());
        }
        dict->name = strdup(enmr.name().c_str());
        dict->metadata = nullptr;
        if (enmr.ordered()) {
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
    arrow_schema->release = &ArrowAdapter::release_schema;
    return arrow_schema;
}

FilterList ArrowAdapter::_create_filter_list(
    std::string filters, std::shared_ptr<Context> ctx) {
    return ArrowAdapter::_create_filter_list(json::parse(filters), ctx);
}

FilterList ArrowAdapter::_create_filter_list(
    json filters, std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    for (auto filter : filters) {
        ArrowAdapter::_append_to_filter_list(filter_list, filter, ctx);
    }

    return filter_list;
}

FilterList ArrowAdapter::_create_attr_filter_list(
    std::string name,
    PlatformConfig platform_config,
    std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    if (platform_config.attrs.empty()) {
        filter_list.add_filter(Filter(*ctx, TILEDB_FILTER_ZSTD));
    } else {
        json attr_options = json::parse(platform_config.attrs);
        if (attr_options.find(name) != attr_options.end() &&
            attr_options[name].find("filters") != attr_options[name].end()) {
            filter_list = ArrowAdapter::_create_filter_list(
                attr_options[name]["filters"], ctx);
        } else {
            filter_list.add_filter(Filter(*ctx, TILEDB_FILTER_ZSTD));
        }
    }

    return filter_list;
}

FilterList ArrowAdapter::_create_dim_filter_list(
    std::string name,
    PlatformConfig platform_config,
    std::string soma_type,
    std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    if (platform_config.dims.empty()) {
        filter_list.add_filter(
            ArrowAdapter::_get_zstd_default(platform_config, soma_type, ctx));
    } else {
        json dim_options = json::parse(platform_config.dims);
        if (dim_options.find(name) != dim_options.end() &&
            dim_options[name].find("filters") != dim_options[name].end()) {
            filter_list = ArrowAdapter::_create_filter_list(
                dim_options[name]["filters"], ctx);
        } else {
            filter_list.add_filter(ArrowAdapter::_get_zstd_default(
                platform_config, soma_type, ctx));
        }
    }

    return filter_list;
}

Filter ArrowAdapter::_get_zstd_default(
    PlatformConfig platform_config,
    std::string soma_type,
    std::shared_ptr<Context> ctx) {
    Filter zstd_filter(*ctx, TILEDB_FILTER_ZSTD);
    if (soma_type == "SOMADataFrame") {
        zstd_filter.set_option(
            TILEDB_COMPRESSION_LEVEL, platform_config.dataframe_dim_zstd_level);
    } else if (soma_type == "SOMASparseNDArray") {
        zstd_filter.set_option(
            TILEDB_COMPRESSION_LEVEL,
            platform_config.sparse_nd_array_dim_zstd_level);
    } else if (soma_type == "SOMADenseNDArray") {
        zstd_filter.set_option(
            TILEDB_COMPRESSION_LEVEL,
            platform_config.dense_nd_array_dim_zstd_level);
    }
    return zstd_filter;
}

void ArrowAdapter::_append_to_filter_list(
    FilterList filter_list, json value, std::shared_ptr<Context> ctx) {
    std::map<std::string, tiledb_filter_type_t> convert_filter = {
        {"GZIP", TILEDB_FILTER_GZIP},
        {"ZSTD", TILEDB_FILTER_ZSTD},
        {"LZ4", TILEDB_FILTER_LZ4},
        {"BZIP2", TILEDB_FILTER_BZIP2},
        {"RLE", TILEDB_FILTER_RLE},
        {"DELTA", TILEDB_FILTER_DELTA},
        {"DOUBLE_DELTA", TILEDB_FILTER_DOUBLE_DELTA},
        {"BIT_WIDTH_REDUCTION", TILEDB_FILTER_BIT_WIDTH_REDUCTION},
        {"BITSHUFFLE", TILEDB_FILTER_BITSHUFFLE},
        {"BYTESHUFFLE", TILEDB_FILTER_BYTESHUFFLE},
        {"POSITIVE_DELTA", TILEDB_FILTER_POSITIVE_DELTA},
        {"CHECKSUM_MD5", TILEDB_FILTER_CHECKSUM_MD5},
        {"CHECKSUM_SHA256", TILEDB_FILTER_CHECKSUM_SHA256},
        {"DICTIONARY_ENCODING", TILEDB_FILTER_DICTIONARY},
        {"SCALE_FLOAT", TILEDB_FILTER_SCALE_FLOAT},
        {"XOR", TILEDB_FILTER_XOR},
        {"WEBP", TILEDB_FILTER_WEBP},
        {"NOOP", TILEDB_FILTER_NONE},
        {"NONE", TILEDB_FILTER_NONE},
    };

    try {
        if (value.is_string()) {
            filter_list.add_filter(Filter(*ctx, convert_filter.at(value)));
        } else {
            Filter filter(*ctx, convert_filter.at(value["name"]));
            for (auto& [key, value] : value.items()) {
                ArrowAdapter::_set_filter_option(filter, key, value);
            }
            filter_list.add_filter(filter);
        }
    } catch (std::out_of_range& e) {
        throw TileDBSOMAError(std::format(
            "Invalid filter {} passed to PlatformConfig", std::string(value)));
    }
}

void ArrowAdapter::_set_filter_option(
    Filter filter, std::string option_name, json value) {
    if (option_name == "name") {
        return;
    }

    std::map<std::string, tiledb_filter_option_t> convert_option = {
        {"COMPRESSION_LEVEL", TILEDB_COMPRESSION_LEVEL},
        {"BIT_WIDTH_MAX_WINDOW", TILEDB_BIT_WIDTH_MAX_WINDOW},
        {"POSITIVE_DELTA_MAX_WINDOW", TILEDB_POSITIVE_DELTA_MAX_WINDOW},
        {"SCALE_FLOAT_BYTEWIDTH", TILEDB_SCALE_FLOAT_BYTEWIDTH},
        {"SCALE_FLOAT_FACTOR", TILEDB_SCALE_FLOAT_FACTOR},
        {"SCALE_FLOAT_OFFSET", TILEDB_SCALE_FLOAT_OFFSET},
        {"WEBP_INPUT_FORMAT", TILEDB_WEBP_INPUT_FORMAT},
        {"WEBP_QUALITY", TILEDB_WEBP_QUALITY},
        {"WEBP_LOSSLESS", TILEDB_WEBP_LOSSLESS},
        {"COMPRESSION_REINTERPRET_DATATYPE",
         TILEDB_COMPRESSION_REINTERPRET_DATATYPE},
    };

    auto option = convert_option[option_name];
    switch (option) {
        case TILEDB_COMPRESSION_LEVEL:
            filter.set_option(option, value.get<int32_t>());
            break;
        case TILEDB_BIT_WIDTH_MAX_WINDOW:
        case TILEDB_POSITIVE_DELTA_MAX_WINDOW:
            filter.set_option(option, value.get<uint32_t>());
            break;
        case TILEDB_SCALE_FLOAT_BYTEWIDTH:
            filter.set_option(option, value.get<uint64_t>());
            break;
        case TILEDB_SCALE_FLOAT_FACTOR:
        case TILEDB_SCALE_FLOAT_OFFSET:
            filter.set_option(option, value.get<double>());
            break;
        case TILEDB_WEBP_QUALITY:
            filter.set_option(option, value.get<float>());
            break;
        case TILEDB_WEBP_INPUT_FORMAT:
        case TILEDB_WEBP_LOSSLESS:
        case TILEDB_COMPRESSION_REINTERPRET_DATATYPE:
            filter.set_option(option, value.get<uint8_t>());
            break;
        default:
            throw TileDBSOMAError(
                std::format("Invalid option {} passed to filter", option_name));
    }
}

Dimension ArrowAdapter::_create_dim(
    tiledb_datatype_t type,
    std::string name,
    const void* buff,
    std::shared_ptr<Context> ctx) {
    switch (type) {
        case TILEDB_STRING_ASCII:
            return Dimension::create(*ctx, name, type, nullptr, nullptr);
        case TILEDB_TIME_SEC:
        case TILEDB_TIME_MS:
        case TILEDB_TIME_US:
        case TILEDB_TIME_NS:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS: {
            // Sadly we cannot put this in the centralized _create_dim_aux
            // in the header file. That's because we need utils/logger.h
            // -- which is a _fixed_ relative path from _this_ .cc file
            // but a _varying_ relative path from all the places that
            // #include arrow_adapter.h. Hence the code duplication in
            // logging statements. :(
            uint64_t* b = (uint64_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return Dimension::create(*ctx, name, type, b, b + 2);
        }
        case TILEDB_INT8: {
            int8_t* b = (int8_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (int8_t*)buff);
        }
        case TILEDB_UINT8: {
            uint8_t* b = (uint8_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (uint8_t*)buff);
        }
        case TILEDB_INT16: {
            int16_t* b = (int16_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (int16_t*)buff);
        }
        case TILEDB_UINT16: {
            uint16_t* b = (uint16_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (uint16_t*)buff);
        }
        case TILEDB_INT32: {
            int32_t* b = (int32_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (int32_t*)buff);
        }
        case TILEDB_UINT32: {
            uint32_t* b = (uint32_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (uint32_t*)buff);
        }
        case TILEDB_INT64: {
            int64_t* b = (int64_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (int64_t*)buff);
        }
        case TILEDB_UINT64: {
            uint64_t* b = (uint64_t*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (uint64_t*)buff);
        }
        case TILEDB_FLOAT32: {
            float* b = (float*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (float*)buff);
        }
        case TILEDB_FLOAT64: {
            double* b = (double*)buff;
            LOG_DEBUG(std::format(
                "_create_dim name={} b={} b1={} b2={}",
                name,
                b[0],
                b[1],
                b[2]));
            return ArrowAdapter::_create_dim_aux(ctx, name, (double*)buff);
        }
        default:
            throw TileDBSOMAError(std::format(
                "ArrowAdapter: Unsupported TileDB dimension: {} ",
                tiledb::impl::type_to_str(type)));
    }
}

tiledb_layout_t ArrowAdapter::_get_order(std::string order) {
    std::transform(
        order.begin(), order.end(), order.begin(), [](unsigned char c) {
            return std::tolower(c);
        });

    std::map<std::string, tiledb_layout_t> convert_order = {
        {"row-major", TILEDB_ROW_MAJOR},
        {"row_major", TILEDB_ROW_MAJOR},
        {"row", TILEDB_ROW_MAJOR},
        {"col-major", TILEDB_COL_MAJOR},
        {"col_major", TILEDB_COL_MAJOR},
        {"column-major", TILEDB_COL_MAJOR},
        {"col", TILEDB_COL_MAJOR},
        {"hilbert", TILEDB_HILBERT},
        {"unordered", TILEDB_UNORDERED},
    };

    try {
        return convert_order[order];
    } catch (const std::out_of_range& e) {
        throw TileDBSOMAError(
            std::format("Invalid order {} passed to PlatformConfig", order));
    }
}

std::tuple<ArraySchema, nlohmann::json>
ArrowAdapter::tiledb_schema_from_arrow_schema(
    std::shared_ptr<Context> ctx,
    const std::unique_ptr<ArrowSchema>& arrow_schema,
    const ArrowTable& index_column_info,
    const std::optional<SOMACoordinateSpace>& coordinate_space,
    std::string soma_type,
    bool is_sparse,
    PlatformConfig platform_config) {
    auto& index_column_array = index_column_info.first;
    auto& index_column_schema = index_column_info.second;

    ArraySchema schema(*ctx, is_sparse ? TILEDB_SPARSE : TILEDB_DENSE);
    Domain domain(*ctx);

    schema.set_capacity(platform_config.capacity);

    if (!platform_config.offsets_filters.empty()) {
        schema.set_offsets_filter_list(ArrowAdapter::_create_filter_list(
            platform_config.offsets_filters, ctx));
    }

    if (!platform_config.validity_filters.empty()) {
        schema.set_validity_filter_list(ArrowAdapter::_create_filter_list(
            platform_config.validity_filters, ctx));
    }

    schema.set_allows_dups(platform_config.allows_duplicates);

    if (platform_config.tile_order) {
        schema.set_tile_order(
            ArrowAdapter::_get_order(*platform_config.tile_order));
    }

    if (platform_config.cell_order) {
        schema.set_cell_order(
            ArrowAdapter::_get_order(*platform_config.cell_order));
    }

    std::vector<std::shared_ptr<SOMAColumn>> columns;

    for (int64_t sch_idx = 0; sch_idx < arrow_schema->n_children; ++sch_idx) {
        auto child = arrow_schema->children[sch_idx];
        std::string_view type_metadata;

        if (ArrowMetadataHasKey(
                child->metadata,
                ArrowCharView(ARROW_DATATYPE_METADATA_KEY.c_str()))) {
            ArrowStringView out;
            NANOARROW_THROW_NOT_OK(ArrowMetadataGetValue(
                child->metadata,
                ArrowCharView(ARROW_DATATYPE_METADATA_KEY.c_str()),
                &out));

            type_metadata = std::string_view(out.data, out.size_bytes);
        }

        LOG_DEBUG(std::format(
            "[ArrowAdapter] schema pass for child {} name '{}'",
            sch_idx,
            std::string(child->name)));

        bool isattr = true;

        for (int64_t i = 0; i < index_column_schema->n_children; ++i) {
            if (strcmp(child->name, index_column_schema->children[i]->name) ==
                0) {
                if (strcmp(child->name, SOMA_GEOMETRY_COLUMN_NAME.c_str()) ==
                    0) {
                    columns.push_back(SOMAGeometryColumn::create(
                        ctx,
                        child,
                        index_column_schema->children[i],
                        index_column_array->children[i],
                        coordinate_space.value(),
                        soma_type,
                        type_metadata,
                        platform_config));
                } else {
                    columns.push_back(SOMADimension::create(
                        ctx,
                        index_column_schema->children[i],
                        index_column_array->children[i],
                        soma_type,
                        type_metadata,
                        platform_config));
                }
                isattr = false;
                LOG_DEBUG(std::format(
                    "[ArrowAdapter] adding dimension {}", child->name));
                break;
            }
        }

        if (isattr) {
            columns.push_back(SOMAAttribute::create(
                ctx, child, type_metadata, platform_config));
            LOG_DEBUG(
                std::format("[ArrowAdapter] adding attribute {}", child->name));
        }
    }

    LOG_DEBUG(std::format("[ArrowAdapter] Additional schema metadata"));
    nlohmann::json soma_schema_extension;
    soma_schema_extension[TILEDB_SOMA_SCHEMA_COL_KEY] = nlohmann::json::array();
    soma_schema_extension["version"] = TILEDB_SOMA_SCHEMA_VERSION;

    // Unit tests expect dimension order should match the index column schema
    // and NOT the Arrow schema
    // We generate the additional schema metadata here to ensure that the
    // serialized column order matches the expected schema order
    for (int64_t i = 0; i < index_column_schema->n_children; ++i) {
        LOG_DEBUG(std::format("[ArrowAdapter] child {}", i));
        const auto column = util::find_column_by_name(
            columns, index_column_schema->children[i]->name);

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
                ArraySchemaExperimental::add_enumeration(
                    *ctx, schema, enumeration);
            }
        }

        if (column->tiledb_attributes().has_value()) {
            auto attributes = column->tiledb_attributes().value();
            for (const auto& attribute : attributes) {
                schema.add_attribute(attribute);
            }
        }
    }

    for (const auto& column : columns | std::views::filter([](const auto& col) {
                                  return !col->isIndexColumn();
                              })) {
        column->serialize(soma_schema_extension[TILEDB_SOMA_SCHEMA_COL_KEY]);

        if (column->tiledb_enumerations().has_value()) {
            auto enumerations = column->tiledb_enumerations().value();
            for (const auto& enumeration : enumerations) {
                ArraySchemaExperimental::add_enumeration(
                    *ctx, schema, enumeration);
            }
        }

        if (column->tiledb_attributes().has_value()) {
            auto attributes = column->tiledb_attributes().value();
            for (const auto& attribute : attributes) {
                schema.add_attribute(attribute);
            }
        }
    }

    LOG_DEBUG(std::format("[ArrowAdapter] set_domain"));
    schema.set_domain(domain);

    LOG_DEBUG(std::format(
        "[ArrowAdapter] index_column_info length {}",
        index_column_array->length));

    // Note: this must be done after we've got the core domain, since the
    // NDRectangle constructor requires access to the core domain.

    CurrentDomain current_domain(*ctx);
    NDRectangle ndrect(*ctx, domain);

    for (auto column : columns) {
        if (!column->isIndexColumn()) {
            continue;
        }

        column->set_current_domain_slot(
            ndrect,
            get_table_any_column_by_name<2>(
                index_column_info, column->name(), 3));

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

    LOG_DEBUG(std::format(
        "[ArrowAdapter] before setting current_domain from ndrect"));
    ArraySchemaExperimental::set_current_domain(*ctx, schema, current_domain);
    LOG_DEBUG(
        std::format("[ArrowAdapter] after setting current_domain from ndrect"));

    LOG_DEBUG(std::format("[ArrowAdapter] check"));
    schema.check();

    LOG_DEBUG(std::format("[ArrowAdapter] returning"));
    return std::make_tuple(schema, soma_schema_extension);
}

Dimension ArrowAdapter::tiledb_dimension_from_arrow_schema(
    std::shared_ptr<Context> ctx,
    ArrowSchema* schema,
    ArrowArray* array,
    std::string soma_type,
    std::string_view type_metadata,
    std::string prefix,
    std::string suffix,
    PlatformConfig platform_config) {
    auto type = ArrowAdapter::to_tiledb_format(schema->format, type_metadata);

    if (ArrowAdapter::arrow_is_var_length_type(schema->format)) {
        type = TILEDB_STRING_ASCII;
    }

    auto col_name = prefix + std::string(schema->name) + suffix;

    FilterList filter_list = ArrowAdapter::_create_dim_filter_list(
        col_name, platform_config, soma_type, ctx);

    if (array->length != 5) {
        throw TileDBSOMAError(std::format(
            "ArrowAdapter: unexpected length {} != 5 for name "
            "'{}'",
            array->length,
            col_name));
    }

    const void* buff = array->buffers[1];
    auto dim = ArrowAdapter::_create_dim(type, col_name, buff, ctx);
    dim.set_filter_list(filter_list);

    return dim;
}

Dimension ArrowAdapter::tiledb_dimension_from_arrow_schema_ext(
    std::shared_ptr<Context> ctx,
    ArrowSchema* schema,
    ArrowArray* array,
    std::string soma_type,
    std::string_view type_metadata,
    std::string prefix,
    std::string suffix,
    PlatformConfig platform_config) {
    if (strcmp(schema->format, "+l") != 0) {
        throw TileDBSOMAError(
            std::format("[tiledb_dimension_from_arrow_schema_ext] Schema "
                        "should be of type list."));
    }

    if (schema->n_children != 1) {
        throw TileDBSOMAError(
            std::format("[tiledb_dimension_from_arrow_schema_ext] Schema "
                        "should have exactly 1 child"));
    }

    auto type = ArrowAdapter::to_tiledb_format(
        schema->children[0]->format, type_metadata);

    if (ArrowAdapter::arrow_is_var_length_type(schema->format)) {
        type = TILEDB_STRING_ASCII;
    }

    auto col_name = prefix + std::string(schema->name) + suffix;

    FilterList filter_list = ArrowAdapter::_create_dim_filter_list(
        col_name, platform_config, soma_type, ctx);

    if (array->length != 5) {
        throw TileDBSOMAError(std::format(
            "ArrowAdapter: unexpected length {} != 5 for name "
            "'{}'",
            array->length,
            col_name));
    }

    const void* buff = array->children[0]->buffers[1];
    auto dim = ArrowAdapter::_create_dim(type, col_name, buff, ctx);
    dim.set_filter_list(filter_list);

    return dim;
}

std::pair<Attribute, std::optional<Enumeration>>
ArrowAdapter::tiledb_attribute_from_arrow_schema(
    std::shared_ptr<Context> ctx,
    ArrowSchema* arrow_schema,
    std::string_view type_metadata,
    PlatformConfig platform_config) {
    auto type = ArrowAdapter::to_tiledb_format(
        arrow_schema->format, type_metadata);

    Attribute attr(*ctx, arrow_schema->name, type);

    FilterList filter_list = ArrowAdapter::_create_attr_filter_list(
        arrow_schema->name, platform_config, ctx);
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
        auto enmr_label = util::get_enmr_label(
            arrow_schema, arrow_schema->dictionary);
        enmr = Enumeration::create_empty(
            *ctx,
            enmr_label,
            enmr_type,
            ArrowAdapter::arrow_is_var_length_type(enmr_format) ?
                TILEDB_VAR_NUM :
                1,
            arrow_schema->flags & ARROW_FLAG_DICTIONARY_ORDERED);
        AttributeExperimental::set_enumeration_name(*ctx, attr, enmr_label);
        LOG_DEBUG(std::format(
            "[ArrowAdapter] dictionary for '{}' as '{}' '{}'",
            std::string(arrow_schema->name),
            tiledb::impl::type_to_str(enmr_type),
            std::string(enmr_format)));
    }

    return {attr, enmr};
}

inline void exitIfError(const ArrowErrorCode ec, const std::string& msg) {
    if (ec != NANOARROW_OK)
        throw TileDBSOMAError(
            std::format("ArrowAdapter: Arrow Error {} ", msg));
}

std::pair<std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>>
ArrowAdapter::to_arrow(std::shared_ptr<ColumnBuffer> column) {
    std::unique_ptr<ArrowSchema> schema = std::make_unique<ArrowSchema>();
    std::unique_ptr<ArrowArray> array = std::make_unique<ArrowArray>();
    auto sch = schema.get();
    auto arr = array.get();

    auto coltype = to_arrow_format(column->type()).data();
    auto natype = to_nanoarrow_type(coltype);
    exitIfError(ArrowSchemaInitFromType(sch, natype), "Bad schema init");
    exitIfError(
        ArrowSchemaSetName(sch, column->name().data()), "Bad schema name");
    exitIfError(
        ArrowSchemaAllocateChildren(sch, 0), "Bad schema children alloc");
    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    schema->release = &release_schema;

    // this will be 3 for char vecs and 2 for enumerations
    int n_buffers = column->is_var() ? 3 : 2;

    // Create an ArrowBuffer to manage the lifetime of `column`.
    // - `arrow_buffer` holds shared_ptr to `column`, increments
    //   the use count and keeps the ColumnBuffer data alive.
    // - When the arrow array is released, `array->release()` is
    //   called with `arrow_buffer` in `private_data`.
    //   `arrow_buffer` is deleted, which decrements the the
    //   `column` use count. When the `column` use count reaches
    //   0, the ColumnBuffer data will be deleted.
    auto arrow_buffer = new ArrowBuffer(column);

    exitIfError(ArrowArrayInitFromType(arr, natype), "Bad array init");
    exitIfError(ArrowArrayAllocateChildren(arr, 0), "Bad array children alloc");
    array->length = column->size();

    LOG_TRACE(std::format(
        "[ArrowAdapter] column type {} name {} nbuf {} {} nullable {}",
        to_arrow_format(column->type()).data(),
        column->name().data(),
        n_buffers,
        array->n_buffers,
        column->is_nullable()));

    if (array->n_buffers != n_buffers) {
        throw TileDBSOMAError(std::format(
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

    LOG_TRACE(std::format(
        "[ArrowAdapter] create array name='{}' use_count={}",
        column->name(),
        column.use_count()));

    array->buffers = (const void**)malloc(sizeof(void*) * n_buffers);
    assert(array->buffers != nullptr);
    array->buffers[0] = nullptr;  // validity addressed below
    array->buffers[n_buffers - 1] = column->data<void*>().data();  // data
    if (n_buffers == 3) {
        array->buffers[1] = column->offsets().data();  // offsets
    }

    if (column->is_nullable()) {
        schema->flags |= ARROW_FLAG_NULLABLE;  // it is also set by default

        // Count nulls
        for (size_t i = 0; i < column->size(); ++i) {
            array->null_count += column->validity()[i] == 0;
        }

        // Convert validity bytemap to a bitmap in place
        column->validity_to_bitmap();
        array->buffers[0] = column->validity().data();
    } else {
        schema->flags &= ~ARROW_FLAG_NULLABLE;  // as ArrowSchemaInitFromType
                                                // leads to NULLABLE set
    }

    if (column->is_ordered()) {
        schema->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
    }

    // Workaround to cast TILEDB_BOOL from uint8 to 1-bit Arrow boolean
    if (column->type() == TILEDB_BOOL) {
        column->data_to_bitmap();
    }

    // Workaround for datetime
    if (column->type() == TILEDB_DATETIME_SEC ||
        column->type() == TILEDB_DATETIME_MS ||
        column->type() == TILEDB_DATETIME_NS) {
        free((void*)schema->format);  // free the 'storage' format
        schema->format = strdup(to_arrow_format(column->type()).data());
    }

    // Workaround for date
    if (column->type() == TILEDB_DATETIME_DAY) {
        free((void*)schema->format);  // free the 'storage' format
        schema->format = strdup(to_arrow_format(column->type()).data());
        // TODO: Put in ColumnBuffer
        size_t n = array->length;
        std::vector<int64_t> indata(n);
        std::memcpy(
            indata.data(), column->data<double>().data(), sizeof(int64_t) * n);
        std::vector<int32_t> vec(n);
        for (size_t i = 0; i < n; i++) {
            vec[i] = static_cast<int32_t>(indata[i]);
        }
        std::memcpy(
            (void*)array->buffers[n_buffers - 1],
            vec.data(),
            sizeof(int32_t) * n);
    }

    auto enmr = column->get_enumeration_info();
    if (enmr.has_value()) {
        auto dict_sch = (ArrowSchema*)malloc(sizeof(ArrowSchema));
        auto dict_arr = (ArrowArray*)malloc(sizeof(ArrowArray));

        auto dcoltype = to_arrow_format(enmr->type(), false).data();
        auto dnatype = to_nanoarrow_type(dcoltype);

        exitIfError(
            ArrowSchemaInitFromType(dict_sch, dnatype), "Bad schema init");
        exitIfError(ArrowSchemaSetName(dict_sch, ""), "Bad schema name");
        exitIfError(
            ArrowSchemaAllocateChildren(dict_sch, 0),
            "Bad schema children alloc");
        dict_sch->release = &release_schema;

        exitIfError(
            ArrowArrayInitFromType(dict_arr, dnatype), "Bad array init");
        exitIfError(
            ArrowArrayAllocateChildren(dict_arr, 0),
            "Bad array children alloc");

        // The release function should be the default nanoarrow release
        // function. The dictionary buffers should be released by the parent
        // array release function that we provide. The arrow array as
        // the owner of new buffers should be responsible for clean-up.
        if (enmr->type() == TILEDB_STRING_ASCII ||
            enmr->type() == TILEDB_STRING_UTF8 || enmr->type() == TILEDB_CHAR ||
            enmr->type() == TILEDB_BLOB) {
            dict_arr->length = _set_var_dictionary_buffers(
                enmr.value(), enmr->context(), dict_arr->buffers);
        } else if (enmr->type() == TILEDB_BOOL) {
            dict_arr->length = _set_bool_dictionary_buffers(
                enmr.value(), enmr->context(), dict_arr->buffers);
        } else {
            dict_arr->length = _set_dictionary_buffers(
                enmr.value(), enmr->context(), dict_arr->buffers);
        }

        schema->dictionary = dict_sch;
        array->dictionary = dict_arr;
    }

    return std::pair(std::move(array), std::move(schema));
}

bool ArrowAdapter::arrow_is_var_length_type(const char* format) {
    return (
        (strcmp(format, "U") == 0) || (strcmp(format, "Z") == 0) ||
        (strcmp(format, "u") == 0) || (strcmp(format, "z") == 0));
}

std::string_view ArrowAdapter::to_arrow_format(
    tiledb_datatype_t tiledb_dtype, bool use_large) {
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
        throw std::out_of_range(std::format(
            "ArrowAdapter: Unsupported TileDB type: {} ",
            tiledb::impl::type_to_str(tiledb_dtype)));
    }
}

tiledb_datatype_t ArrowAdapter::to_tiledb_format(
    std::string_view arrow_dtype, std::string_view arrow_dtype_metadata) {
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
        } else if (
            dtype == TILEDB_STRING_UTF8 &&
            arrow_dtype_metadata.compare("WKT") == 0) {
            dtype = TILEDB_GEOM_WKT;
        }

        return dtype;
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(std::format(
            "ArrowAdapter: Unsupported Arrow type: {} ", arrow_dtype));
    }
}

// FIXME: Add more types, maybe make it a map
enum ArrowType ArrowAdapter::to_nanoarrow_type(std::string_view sv) {
    if (sv == "i")
        return NANOARROW_TYPE_INT32;
    else if (sv == "c")
        return NANOARROW_TYPE_INT8;
    else if (sv == "C")
        return NANOARROW_TYPE_UINT8;
    else if (sv == "s")
        return NANOARROW_TYPE_INT16;
    else if (sv == "S")
        return NANOARROW_TYPE_UINT16;
    else if (sv == "I")
        return NANOARROW_TYPE_UINT32;
    else if (sv == "l")
        return NANOARROW_TYPE_INT64;
    else if (sv == "L")
        return NANOARROW_TYPE_UINT64;
    else if (sv == "f")
        return NANOARROW_TYPE_FLOAT;
    else if (sv == "g")
        return NANOARROW_TYPE_DOUBLE;
    else if (sv == "u")
        return NANOARROW_TYPE_STRING;
    else if (sv == "U")
        return NANOARROW_TYPE_LARGE_STRING;
    else if (sv == "b")
        return NANOARROW_TYPE_BOOL;
    else if (sv == "tss:")
        return NANOARROW_TYPE_INT64;  // NB time resolution set indepedently
    else if (sv == "tsm:")
        return NANOARROW_TYPE_INT64;  // NB time resolution set indepedently
    else if (sv == "tsn:")
        return NANOARROW_TYPE_INT64;  // NB time resolution set indepedently
    else if (sv == "tsu:")
        return NANOARROW_TYPE_INT64;  // NB time resolution set indepedently
    else if (sv == "tdD")
        return NANOARROW_TYPE_INT32;  // R Date: fractional days since epoch
    else if (sv == "z")
        return NANOARROW_TYPE_BINARY;
    else if (sv == "Z")
        return NANOARROW_TYPE_LARGE_BINARY;
    else
        throw TileDBSOMAError(
            std::format("ArrowAdapter: Unsupported Arrow format: {} ", sv));
}

std::unique_ptr<ArrowSchema> ArrowAdapter::make_arrow_schema(
    const std::vector<std::string>& names,
    const std::vector<tiledb_datatype_t>& tiledb_datatypes) {
    auto num_names = names.size();
    auto num_types = tiledb_datatypes.size();

    if (num_names != num_types) {
        throw TileDBSOMAError(std::format(
            "ArrowAdapter::make_arrow_schema: internal coding error: num_types "
            "{} != num_names {}",
            num_names,
            num_types));
    }

    auto arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = "+s";  // structure, i.e. non-leaf node
    arrow_schema->name = strdup("parent");
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    arrow_schema->n_children = num_names;  // non-leaf node
    arrow_schema->children = (ArrowSchema**)malloc(
        arrow_schema->n_children * sizeof(ArrowSchema*));
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;

    LOG_DEBUG(std::format(
        "[ArrowAdapter] make_arrow_schema n_children {}",
        arrow_schema->n_children));

    for (int i = 0; i < (int)num_names; i++) {
        ArrowSchema* dim_schema = (ArrowSchema*)malloc(sizeof(ArrowSchema));
        auto arrow_type_name = ArrowAdapter::tdb_to_arrow_type(
            tiledb_datatypes[i]);
        dim_schema->name = strdup(names[i].c_str());
        dim_schema->format = strdup(arrow_type_name.c_str());
        dim_schema->metadata = nullptr;
        dim_schema->flags = 0;
        dim_schema->n_children = 0;      // leaf node
        dim_schema->children = nullptr;  // leaf node
        dim_schema->dictionary = nullptr;
        dim_schema->release = &ArrowAdapter::release_schema;
        dim_schema->private_data = nullptr;

        arrow_schema->children[i] = dim_schema;
        LOG_TRACE(std::format(
            "[ArrowAdapter] make_arrow_schema child {} format {} name {}",
            i,
            dim_schema->format,
            dim_schema->name));

        if (strcmp(dim_schema->name, SOMA_GEOMETRY_COLUMN_NAME.c_str()) == 0) {
            nanoarrow::UniqueBuffer buffer;
            ArrowMetadataBuilderInit(buffer.get(), nullptr);
            ArrowMetadataBuilderAppend(
                buffer.get(),
                ArrowCharView(ARROW_DATATYPE_METADATA_KEY.c_str()),
                ArrowCharView(
                    tiledb_datatypes[i] == TILEDB_GEOM_WKB ? "WKB" : "WKT"));
            ArrowSchemaSetMetadata(
                dim_schema,
                std::string((char*)buffer->data, buffer->size_bytes).c_str());
        }
    }

    return arrow_schema;
}

std::unique_ptr<ArrowSchema> ArrowAdapter::make_arrow_schema_parent(
    size_t num_columns, std::string_view name) {
    auto arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = strdup("+s");  // structure, i.e. non-leaf node
    arrow_schema->name = strdup(name.data());
    arrow_schema->metadata = nullptr;
    arrow_schema->flags = 0;
    arrow_schema->n_children = static_cast<int64_t>(
        num_columns);  // non-leaf node
    arrow_schema->children = (ArrowSchema**)malloc(
        arrow_schema->n_children * sizeof(ArrowSchema*));
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->private_data = nullptr;

    for (size_t i = 0; i < num_columns; i++) {
        arrow_schema->children[i] = nullptr;
    }

    LOG_DEBUG(std::format(
        "[ArrowAdapter] make_arrow_schema n_children {}",
        arrow_schema->n_children));

    return arrow_schema;
}

std::unique_ptr<ArrowArray> ArrowAdapter::make_arrow_array_parent(
    size_t num_columns) {
    auto arrow_array = std::make_unique<ArrowArray>();

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

    arrow_array->children = (ArrowArray**)malloc(
        num_columns * sizeof(ArrowArray*));
    for (size_t i = 0; i < num_columns; i++) {
        arrow_array->children[i] = nullptr;
    }

    LOG_DEBUG(std::format(
        "[ArrowAdapter] make_arrow_array n_children {}",
        arrow_array->n_children));

    return arrow_array;
}

void ArrowAdapter::log_make_arrow_array_child(ArrowArray* child) {
    LOG_TRACE(std::format(
        "[ArrowAdapter] make_arrow_array_child length {} n_buffers {}",
        child->length,
        child->n_buffers));
}

void ArrowAdapter::_check_shapes(
    ArrowArray* arrow_array, ArrowSchema* arrow_schema) {
    if (arrow_array->n_children != arrow_schema->n_children) {
        throw std::runtime_error(
            "ArrowAdapter::_check_shapes: internal coding error: data/schema "
            "mismatch");
    }
    for (int64_t i = 0; i < arrow_array->n_children; i++) {
        _check_shapes(arrow_array->children[i], arrow_schema->children[i]);
    }
}

int64_t ArrowAdapter::_get_column_index_from_name(
    const ArrowTable& arrow_table, std::string column_name) {
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

    throw std::runtime_error(std::format(
        "ArrowAdapter::_check_shapes: column {} not found", column_name));
}

ArrowArray* ArrowAdapter::_get_and_check_column(
    const ArrowTable& arrow_table,
    int64_t column_index,
    int64_t expected_n_buffers) {
    ArrowArray* arrow_array = arrow_table.first.get();
    if (column_index < 0 || column_index >= arrow_array->n_children) {
        throw std::runtime_error(std::format(
            "ArrowAdapter::_get_and_check_column: column index {} out of "
            "bounds {}..{}",
            column_index,
            0,
            arrow_array->n_children - 1));
    }

    ArrowArray* child = arrow_array->children[column_index];

    if (child->n_children != 0) {
        throw std::runtime_error(std::format(
            "ArrowAdapter::_get_and_check_column: column index {} is "
            "non-terminal",
            column_index));
    }

    if (expected_n_buffers == 2) {
        if (child->n_buffers != 2) {
            throw std::runtime_error(std::format(
                "ArrowAdapter::_get_and_check_column: column index {} "
                "has buffer count {}; expected 2 for non-string data",
                column_index,
                child->n_buffers));
        }

    } else if (expected_n_buffers == 3) {
        if (child->n_buffers != 3) {
            throw std::runtime_error(std::format(
                "ArrowAdapter::_get_and_check_column: column index {} is "
                "has buffer count {}; expected 3 for string data",
                column_index,
                child->n_buffers));
        }

    } else {
        throw std::runtime_error(std::format(
            "ArrowAdapter::_get_and_check_column: internal coding error: "
            "expected_n_buffers {} is "
            "neither 2 nor 3.",
            expected_n_buffers));
    }

    return child;
}

std::unique_ptr<ArrowArray> ArrowAdapter::arrow_array_insert_at_index(
    std::unique_ptr<ArrowArray> parent_array,
    std::vector<std::unique_ptr<ArrowArray>> child_arrays,
    int64_t index) {
    if (parent_array->n_children < index || index < 0) {
        throw std::runtime_error(
            "[ArrowAdapter][arrow_array_insert_at_index] Invalid index to "
            "insert array");
    }

    if (child_arrays.size() == 0) {
        return parent_array;
    }

    auto array = make_arrow_array_parent(
        static_cast<size_t>(parent_array->n_children) + child_arrays.size());

    for (int64_t i = 0; i < array->n_children; ++i) {
        int64_t idx = i <= index ? i : i - child_arrays.size();
        array->children[i] = (ArrowArray*)malloc(sizeof(ArrowArray));

        if (i >= index &&
            i < index + static_cast<int64_t>(child_arrays.size())) {
            ArrowArrayMove(
                child_arrays[i - index].release(), array->children[i]);
        } else {
            ArrowArrayMove(parent_array->children[idx], array->children[i]);
        }
    }

    parent_array->release(parent_array.get());

    return array;
}

std::unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_insert_at_index(
    std::unique_ptr<ArrowSchema> parent_schema,
    std::vector<std::unique_ptr<ArrowSchema>> child_schemas,
    int64_t index) {
    if (parent_schema->n_children < index || index < 0) {
        throw std::runtime_error(
            "[ArrowAdapter][arrow_schema_insert_at_index] Invalid index to "
            "insert schema");
    }

    if (child_schemas.size() == 0) {
        return parent_schema;
    }

    auto schema = make_arrow_schema_parent(
        parent_schema->n_children + child_schemas.size());

    for (int64_t i = 0; i < schema->n_children; ++i) {
        int64_t idx = i <= index ? i : i - child_schemas.size();
        schema->children[i] = (ArrowSchema*)malloc(sizeof(ArrowSchema));

        if (i >= index &&
            i < index + static_cast<int64_t>(child_schemas.size())) {
            ArrowSchemaMove(
                child_schemas[i - index].release(), schema->children[i]);
        } else {
            ArrowSchemaMove(parent_schema->children[idx], schema->children[i]);
        }
    }

    parent_schema->release(parent_schema.get());

    return schema;
}

std::unique_ptr<ArrowArray> ArrowAdapter::arrow_array_remove_at_index(
    std::unique_ptr<ArrowArray> array, int64_t index) {
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

    array->release(array.get());

    return array_new;
}

std::unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_remove_at_index(
    std::unique_ptr<ArrowSchema> schema, int64_t index) {
    if (schema->n_children <= index || index < 0) {
        throw std::runtime_error(
            "[ArrowAdapter][arrow_schema_remove_at_index] Invalid index to "
            "remove child schema");
    }

    auto schema_new = make_arrow_schema_parent(schema->n_children - 1);

    for (int64_t i = 0; i < schema->n_children; ++i) {
        int64_t idx = i <= index ? i : i - 1;

        if (i != index) {
            schema_new->children[idx] = (ArrowSchema*)malloc(
                sizeof(ArrowSchema));
            ArrowSchemaMove(schema->children[i], schema_new->children[idx]);
        }
    }

    schema->release(schema.get());

    return schema_new;
}

size_t ArrowAdapter::_set_var_dictionary_buffers(
    Enumeration& enumeration, const Context& ctx, const void** buffers) {
    const void* data;
    uint64_t data_size;

    ctx.handle_error(tiledb_enumeration_get_data(
        ctx.ptr().get(), enumeration.ptr().get(), &data, &data_size));

    const void* offsets;
    uint64_t offsets_size;
    ctx.handle_error(tiledb_enumeration_get_offsets(
        ctx.ptr().get(), enumeration.ptr().get(), &offsets, &offsets_size));

    size_t count = offsets_size / sizeof(uint64_t);

    std::span<const uint64_t> offsets_v(
        static_cast<const uint64_t*>(offsets), count);

    uint32_t* small_offsets = static_cast<uint32_t*>(
        malloc((count + 1) * sizeof(uint32_t)));
    buffers[2] = malloc(data_size);

    std::memcpy(const_cast<void*>(buffers[2]), data, data_size);
    for (size_t i = 0; i < count; ++i) {
        small_offsets[i] = static_cast<uint32_t>(offsets_v[i]);
    }
    small_offsets[count] = static_cast<uint32_t>(data_size);
    buffers[1] = small_offsets;

    return count;
}

size_t ArrowAdapter::_set_dictionary_buffers(
    Enumeration& enumeration, const Context& ctx, const void** buffers) {
    const void* data;
    uint64_t data_size;

    ctx.handle_error(tiledb_enumeration_get_data(
        ctx.ptr().get(), enumeration.ptr().get(), &data, &data_size));

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
        case TILEDB_INT64:
            return data_size / sizeof(int64_t);
        case TILEDB_UINT64:
            return data_size / sizeof(uint64_t);
        case TILEDB_FLOAT32:
            return data_size / sizeof(float_t);
        case TILEDB_FLOAT64:
            return data_size / sizeof(double_t);
        default:
            throw TileDBSOMAError(std::format(
                "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                tiledb::impl::type_to_str(enumeration.type())));
    }
}

size_t ArrowAdapter::_set_bool_dictionary_buffers(
    Enumeration& enumeration, const Context& ctx, const void** buffers) {
    const void* data;
    uint64_t data_size;

    ctx.handle_error(tiledb_enumeration_get_data(
        ctx.ptr().get(), enumeration.ptr().get(), &data, &data_size));

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
