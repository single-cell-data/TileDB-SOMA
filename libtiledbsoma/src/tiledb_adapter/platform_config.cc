/**
 * @file   platform_config.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines classes and helper methods for the TileDB platform configuration.
 */

#include "platform_config.h"
#include "nlohmann/json.hpp"

#include "../utils/common.h"
#include "common/logging/impl/logger.h"

using json = nlohmann::json;

namespace tiledbsoma::utils {

tiledb_layout_t get_order_from_string(std::string order) {
    std::transform(order.begin(), order.end(), order.begin(), [](unsigned char c) { return std::tolower(c); });

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
        throw TileDBSOMAError(fmt::format("Invalid order {} passed to PlatformConfig", order));
    }
}

ArraySchema create_base_tiledb_schema(
    std::shared_ptr<Context> ctx,
    const PlatformConfig& platform_config,
    bool is_sparse,
    std::optional<std::pair<int64_t, int64_t>> timestamp_range) {
    tiledb_array_schema_t* c_schema;
    if (timestamp_range.has_value() && timestamp_range.value().first != 0) {
        ctx->handle_error(tiledb_array_schema_alloc_at_timestamp(
            ctx->ptr().get(), is_sparse ? TILEDB_SPARSE : TILEDB_DENSE, timestamp_range.value().first, &c_schema));
    } else {
        ctx->handle_error(
            tiledb_array_schema_alloc(ctx->ptr().get(), is_sparse ? TILEDB_SPARSE : TILEDB_DENSE, &c_schema));
    }
    ArraySchema schema(*ctx, c_schema);

    schema.set_capacity(platform_config.capacity);

    if (!platform_config.offsets_filters.empty()) {
        schema.set_offsets_filter_list(utils::create_filter_list(platform_config.offsets_filters, ctx));
    }

    if (!platform_config.validity_filters.empty()) {
        schema.set_validity_filter_list(utils::create_filter_list(platform_config.validity_filters, ctx));
    }

    schema.set_allows_dups(platform_config.allows_duplicates);

    if (platform_config.tile_order) {
        schema.set_tile_order(get_order_from_string(*platform_config.tile_order));
    }

    if (platform_config.cell_order) {
        schema.set_cell_order(get_order_from_string(*platform_config.cell_order));
    }

    return schema;
}

json get_filter_list_json(tiledb::FilterList filter_list) {
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
        {TILEDB_COMPRESSION_REINTERPRET_DATATYPE, "COMPRESSION_REINTERPRET_DATATYPE"},
    };

    json filter_list_as_json = {};
    for (uint32_t i = 0; i < filter_list.nfilters(); ++i) {
        json filter_as_json = {};

        auto filter = filter_list.filter(i);
        filter_as_json.emplace("name", tiledb::Filter::to_str(filter.filter_type()));

        switch (filter.filter_type()) {
            case TILEDB_FILTER_GZIP:
            case TILEDB_FILTER_ZSTD:
            case TILEDB_FILTER_LZ4:
            case TILEDB_FILTER_BZIP2:
            case TILEDB_FILTER_RLE:
            case TILEDB_FILTER_DICTIONARY:
                filter_as_json.emplace("COMPRESSION_LEVEL", filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL));
                break;

            case TILEDB_FILTER_DELTA:
            case TILEDB_FILTER_DOUBLE_DELTA:
                filter_as_json.emplace("COMPRESSION_LEVEL", filter.get_option<int32_t>(TILEDB_COMPRESSION_LEVEL));
                filter_as_json.emplace(
                    "COMPRESSION_REINTERPRET_DATATYPE",
                    filter.get_option<uint8_t>(TILEDB_COMPRESSION_REINTERPRET_DATATYPE));
                break;

            case TILEDB_FILTER_BIT_WIDTH_REDUCTION:
                filter_as_json.emplace(
                    "BIT_WIDTH_MAX_WINDOW", filter.get_option<uint32_t>(TILEDB_BIT_WIDTH_MAX_WINDOW));
                break;

            case TILEDB_FILTER_POSITIVE_DELTA:
                filter_as_json.emplace(
                    "POSITIVE_DELTA_MAX_WINDOW", filter.get_option<uint32_t>(TILEDB_POSITIVE_DELTA_MAX_WINDOW));
                break;

            case TILEDB_FILTER_SCALE_FLOAT:
                filter_as_json.emplace("SCALE_FLOAT_FACTOR", filter.get_option<double>(TILEDB_SCALE_FLOAT_FACTOR));
                filter_as_json.emplace("SCALE_FLOAT_OFFSET", filter.get_option<double>(TILEDB_SCALE_FLOAT_OFFSET));
                filter_as_json.emplace(
                    "SCALE_FLOAT_BYTEWIDTH", filter.get_option<uint64_t>(TILEDB_SCALE_FLOAT_BYTEWIDTH));
                break;

            case TILEDB_FILTER_WEBP:
                filter_as_json.emplace("WEBP_INPUT_FORMAT", filter.get_option<uint8_t>(TILEDB_WEBP_INPUT_FORMAT));
                filter_as_json.emplace("WEBP_QUALITY", filter.get_option<float>(TILEDB_WEBP_QUALITY));
                filter_as_json.emplace("WEBP_LOSSLESS", filter.get_option<uint8_t>(TILEDB_WEBP_LOSSLESS));
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
            default:
                throw TileDBSOMAError(
                    fmt::format(
                        "Internal error: unrecognized filter type '{}'", tiledb::Filter::to_str(filter.filter_type())));
        }
        filter_list_as_json.emplace_back(filter_as_json);
    }
    return filter_list_as_json;
}

json get_attrs_filter_list_json(const tiledb::ArraySchema& tiledb_schema) {
    json attrs_filter_list_as_json;
    for (const auto& attr : tiledb_schema.attributes()) {
        json attr_info = {{"filters", get_filter_list_json(attr.second.filter_list())}};
        attrs_filter_list_as_json.emplace(attr.first, attr_info);
    }
    return attrs_filter_list_as_json;
}

json get_dims_list_json(const tiledb::ArraySchema& tiledb_schema) {
    json dims_as_json;
    for (const auto& dim : tiledb_schema.domain().dimensions()) {
        json dim_info = {{"tile", dim.tile_extent_to_str()}, {"filters", get_filter_list_json(dim.filter_list())}};
        dims_as_json.emplace(dim.name(), dim_info);
    }
    return dims_as_json;
}

void set_filter_option(Filter filter, std::string option_name, json value) {
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
        {"COMPRESSION_REINTERPRET_DATATYPE", TILEDB_COMPRESSION_REINTERPRET_DATATYPE},
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
            throw TileDBSOMAError(fmt::format("Invalid option {} passed to filter", option_name));
    }
}

void append_to_filter_list(FilterList filter_list, json value, std::shared_ptr<Context> ctx) {
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
                set_filter_option(filter, key, value);
            }
            filter_list.add_filter(filter);
        }
    } catch (std::out_of_range& e) {
        throw TileDBSOMAError(fmt::format("Invalid filter {} passed to PlatformConfig", std::string(value)));
    }
}

FilterList create_filter_list(json filters, std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    for (auto filter : filters) {
        append_to_filter_list(filter_list, filter, ctx);
    }

    return filter_list;
}

FilterList create_filter_list(std::string filters, std::shared_ptr<Context> ctx) {
    return create_filter_list(json::parse(filters), ctx);
}

FilterList create_attr_filter_list(std::string name, PlatformConfig platform_config, std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    if (platform_config.attrs.empty()) {
        filter_list.add_filter(Filter(*ctx, TILEDB_FILTER_ZSTD));
    } else {
        json attr_options = json::parse(platform_config.attrs);
        if (attr_options.find(name) != attr_options.end() &&
            attr_options[name].find("filters") != attr_options[name].end()) {
            filter_list = create_filter_list(attr_options[name]["filters"], ctx);
        } else {
            filter_list.add_filter(Filter(*ctx, TILEDB_FILTER_ZSTD));
        }
    }

    return filter_list;
}

Filter get_zstd_default(PlatformConfig platform_config, std::string soma_type, std::shared_ptr<Context> ctx) {
    Filter zstd_filter(*ctx, TILEDB_FILTER_ZSTD);
    if (soma_type == "SOMADataFrame") {
        zstd_filter.set_option(TILEDB_COMPRESSION_LEVEL, platform_config.dataframe_dim_zstd_level);
    } else if (soma_type == "SOMASparseNDArray") {
        zstd_filter.set_option(TILEDB_COMPRESSION_LEVEL, platform_config.sparse_nd_array_dim_zstd_level);
    } else if (soma_type == "SOMADenseNDArray") {
        zstd_filter.set_option(TILEDB_COMPRESSION_LEVEL, platform_config.dense_nd_array_dim_zstd_level);
    }
    return zstd_filter;
}

FilterList create_dim_filter_list(
    std::string name, PlatformConfig platform_config, std::string soma_type, std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    if (platform_config.dims.empty()) {
        filter_list.add_filter(get_zstd_default(platform_config, soma_type, ctx));
    } else {
        json dim_options = json::parse(platform_config.dims);
        if (dim_options.find(name) != dim_options.end() &&
            dim_options[name].find("filters") != dim_options[name].end()) {
            filter_list = create_filter_list(dim_options[name]["filters"], ctx);
        } else {
            filter_list.add_filter(get_zstd_default(platform_config, soma_type, ctx));
        }
    }

    return filter_list;
}

PlatformConfig platform_config_from_tiledb_schema(tiledb::ArraySchema tiledb_schema) {
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
    platform_config.offsets_filters = get_filter_list_json(tiledb_schema.offsets_filter_list()).dump();
    platform_config.validity_filters = get_filter_list_json(tiledb_schema.validity_filter_list()).dump();
    platform_config.attrs = get_attrs_filter_list_json(tiledb_schema).dump();
    platform_config.dims = get_dims_list_json(tiledb_schema).dump();

    return platform_config;
}

PlatformSchemaConfig platform_schema_config_from_tiledb(tiledb::ArraySchema tiledb_schema) {
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
    platform_config.offsets_filters = get_filter_list_json(tiledb_schema.offsets_filter_list()).dump();
    platform_config.validity_filters = get_filter_list_json(tiledb_schema.validity_filter_list()).dump();
    platform_config.attrs = get_attrs_filter_list_json(tiledb_schema).dump();
    platform_config.dims = get_dims_list_json(tiledb_schema).dump();

    return platform_config;
}
}  // namespace tiledbsoma::utils
