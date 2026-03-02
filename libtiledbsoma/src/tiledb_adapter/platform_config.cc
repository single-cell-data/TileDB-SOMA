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

#include "../utils/common.h"
#include "../utils/util.h"
#include "common/arrow/utils.h"
#include "common/logging/impl/logger.h"

#include <algorithm>
#include <limits>
#include <ranges>
#include <tiledb/tiledb_experimental>
#include <unordered_map>

#include "../soma/soma_attribute.h"
#include "../soma/soma_dimension.h"

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

ArraySchema create_nd_array_schema(
    std::string_view soma_type,
    bool is_sparse,
    std::string_view format,
    std::span<const int64_t> shape,
    std::shared_ptr<tiledb::Context> ctx,
    PlatformConfig platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    tiledb::ArraySchema schema = utils::create_base_tiledb_schema(ctx, platform_config, is_sparse, timestamp);
    tiledb::Domain domain(*ctx);
    int64_t default_extent = is_sparse ? 2048 : std::max(2048 >> shape.size(), 4);

    for (size_t i = 0; i < shape.size(); ++i) {
        if (shape[i] <= 0) {
            throw std::range_error("[create_nd_array_schema] Shape slots must be at least 1");
        }

        std::string name = fmt::format("soma_dim_{}", i);
        int64_t extent = utils::get_dim_extent<int64_t>(
            name, platform_config, default_extent, std::numeric_limits<int64_t>::max());

        tiledb::Dimension dimension = tiledb::Dimension::create<int64_t>(
            *ctx, name, {{0, std::numeric_limits<int64_t>::max() - extent - 1}}, extent);
        dimension.set_filter_list(utils::create_dim_filter_list(name, platform_config, soma_type.data(), ctx));
        domain.add_dimension(dimension);
    }

    schema.set_domain(domain);

    tiledb::Attribute attribute = tiledb::Attribute::create(*ctx, "soma_data", common::arrow::to_tiledb_format(format));
    attribute.set_filter_list(utils::create_attr_filter_list("soma_data", platform_config, ctx));
    schema.add_attribute(attribute);

    tiledb::CurrentDomain current_domain(*ctx);
    tiledb::NDRectangle rect(*ctx, domain);

    for (size_t i = 0; i < shape.size(); ++i) {
        rect.set_range<int64_t>(i, 0, shape[i] - 1);
    }

    current_domain.set_ndrectangle(rect);
    tiledb::ArraySchemaExperimental::set_current_domain(*ctx, schema, current_domain);
    schema.check();

    return schema;
}

ArraySchema create_dataframe_schema(
    std::string_view soma_type,
    ArrowSchema* arrow_schema,
    std::span<const std::string> index_column_names,
    std::span<const DomainRange> index_column_domains,
    std::shared_ptr<tiledb::Context> ctx,
    PlatformConfig platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    tiledb::ArraySchema schema = utils::create_base_tiledb_schema(ctx, platform_config, true, timestamp);
    tiledb::Domain domain(*ctx);

    if (index_column_names.empty()) {
        throw std::range_error("At least one non-index column must be provided");
    }

    if (arrow_schema->n_children) {
        for (const auto& name : index_column_names) {
            auto column_occurrences = std::count_if(
                arrow_schema->children,
                arrow_schema->children + arrow_schema->n_children,
                [&](const auto child_schema) { return name == child_schema->name; });

            if (column_occurrences == 0) {
                // 'soma_joinid` will be added in the schema by default if missing
                if (name != SOMA_JOINID) {
                    throw std::range_error(fmt::format("Missing index column '{}' from schema", name));
                }
            } else if (column_occurrences != 1) {
                throw std::range_error(fmt::format("Multiple columns found in schema for index column '{}'", name));
            }
        }
    }

    if (index_column_names.size() != index_column_domains.size()) {
        throw std::range_error(
            fmt::format(
                "Size mismatch between list of index column names and domains; {} != {}",
                index_column_names.size(),
                index_column_domains.size()));
    }

    std::unordered_map<std::string, DomainRange> index_columns;
    for (size_t i = 0; i < index_column_names.size(); ++i) {
        if (index_columns.contains(index_column_names[i])) {
            throw std::range_error(fmt::format("Duplicate index column found '{}'", index_column_names[i]));
        }

        index_columns.emplace(index_column_names[i], index_column_domains[i]);
    }

    std::vector<std::shared_ptr<SOMAColumn>> columns;

    for (int64_t i = 0; i < arrow_schema->n_children; ++i) {
        std::string_view name(arrow_schema->children[i]->name);
        std::string_view format(arrow_schema->children[i]->format);

        // Verify no illegal use of soma_ prefix
        if (!platform_config.override_naming_restriction && name.starts_with("soma_")) {
            if (name != SOMA_JOINID) {
                throw std::range_error(
                    fmt::format("DataFrame schema may not contain fields with name prefix 'soma_': got '{}'", name));
            }

            if (name == SOMA_JOINID) {
                if (common::arrow::to_tiledb_format(format) != TILEDB_INT64) {
                    throw std::range_error(
                        fmt::format(
                            "'{}' field must be of type Arrow int64 but is {}",
                            SOMA_JOINID,
                            common::arrow::to_arrow_readable(format)));
                }
            }
        }

        if (index_columns.contains(name.data())) {
            if (arrow_schema->dictionary != nullptr) {
                std::range_error(
                    fmt::format(
                        "Cannot set index column '{}' to an enumeration. Index columns do not support enumerations.",
                        name));
            }

            if (name == SOMA_JOINID) {
                columns.push_back(SOMADimension::create_soma_joinid(ctx, "SOMADataFrame", platform_config));
            } else {
                columns.push_back(
                    SOMADimension::create(ctx, arrow_schema->children[i], "SOMADataFrame", platform_config));
            }
        } else {
            columns.push_back(SOMAAttribute::create(ctx, arrow_schema->children[i], platform_config));
        }
    }

    // Inject missing required columns
    if (std::none_of(
            columns.cbegin(), columns.cend(), [](const auto& column) { return column->name() == SOMA_JOINID; })) {
        if (index_columns.contains(SOMA_JOINID.data())) {
            columns.push_back(SOMADimension::create_soma_joinid(ctx, "SOMADataFrame", platform_config));
        } else {
            columns.push_back(
                std::make_shared<SOMAAttribute>(tiledb::Attribute::create<int64_t>(
                    *ctx,
                    SOMA_JOINID.data(),
                    utils::create_attr_filter_list(SOMA_JOINID.data(), platform_config, ctx))));
        }
    }

    // Unit tests expect dimension order should match the index column schema
    // and NOT the Arrow schema
    // We generate the additional schema metadata here to ensure that the
    // serialized column order matches the expected schema order
    for (const auto& column_name : index_column_names) {
        // If a column is specified as `index column` but it isn't present in the schema the following call is expected to throw
        const auto column = util::find_column_by_name(columns, column_name);

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

    schema.set_domain(domain);

    tiledb::CurrentDomain current_domain(*ctx);
    tiledb::NDRectangle rect(*ctx, domain);

    auto decode_domain = []<typename T>(
                             std::shared_ptr<SOMAColumn> column, std::optional<std::pair<T, T>> domain) -> std::any {
        if constexpr (std::is_same_v<T, std::string>) {
            auto current_domain = domain.value_or(std::make_pair<T, T>("", ""));

            if (current_domain.first != "" || current_domain.second != "") {
                throw std::range_error("TileDB str and bytes index-column types do not support domain specification");
            }

            return std::make_any<std::array<std::string, 2>>(std::array<std::string, 2>{{"", ""}});
        } else {
            auto current_domain = domain.value_or(std::make_pair<T, T>(0, 0));

            if (column->name() == SOMA_JOINID) {
                if (current_domain.first < 0) {
                    throw std::range_error(
                        fmt::format(
                            "'{}' indices cannot be negative; got lower bound {}", SOMA_JOINID, current_domain.first));
                } else if (current_domain.second < 0) {
                    throw std::range_error(
                        fmt::format(
                            "'{}' indices cannot be negative; got bound bound {}", SOMA_JOINID, current_domain.second));
                }
            }

            return std::make_any<std::array<T, 2>>(std::array<T, 2>{{current_domain.first, current_domain.second}});
        }
    };

    for (const auto& column : columns | std::views::filter([](const auto& col) { return col->isIndexColumn(); })) {
        std::visit(
            [&](auto&& domain) {
                using T = std::decay_t<decltype(domain)>;

                if constexpr (std::is_same_v<T, std::monostate>) {
                    throw std::runtime_error(
                        fmt::format(
                            "Index column '{}' has no domain attached. This can be caused by the column not being "
                            "present in the schema provided.",
                            column->name()));
                } else {
                    using E = std::decay_t<decltype(domain)>::value_type::first_type;

                    column->set_current_domain_slot(rect, {{decode_domain.template operator()<E>(column, domain)}});
                }
            },
            index_columns[column->name()]);
    }

    current_domain.set_ndrectangle(rect);
    tiledb::ArraySchemaExperimental::set_current_domain(*ctx, schema, current_domain);
    schema.check();

    return schema;
}
}  // namespace tiledbsoma::utils
