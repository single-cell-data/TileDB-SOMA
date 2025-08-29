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

using json = nlohmann::json;

namespace tiledbsoma::utils {

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
