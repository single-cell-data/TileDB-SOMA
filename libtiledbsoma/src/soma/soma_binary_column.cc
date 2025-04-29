/**
 * @file   soma_binary_column.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMABinaryColumn class.
 */

#include "soma_binary_column.h"
#include "../utils/logger.h"

namespace tiledbsoma {
std::shared_ptr<SOMAColumn> SOMABinaryColumn::deserialize(
    const nlohmann::json& soma_schema,
    const Context& ctx,
    const Array& array,
    const std::map<std::string, tiledbsoma::MetadataValue>&) {
    if (soma_schema.contains(TILEDB_SOMA_SCHEMA_COL_ATTR_KEY)) {
        std::vector<std::string>
            attribute_names = soma_schema[TILEDB_SOMA_SCHEMA_COL_ATTR_KEY]
                                  .template get<std::vector<std::string>>();

        if (attribute_names.size() != 1) {
            throw TileDBSOMAError(fmt::format(
                "[SOMABinaryColumn][deserialize] Invalid number of attributes. "
                "Epected 1, got {}",
                attribute_names.size()));
        }

        if (!array.schema().has_attribute(attribute_names[0])) {
            // Attribute probably dropped so skip column reconstruction.
            return nullptr;
        }

        auto attribute = array.schema().attribute(attribute_names[0]);
        auto enumeration_name = AttributeExperimental::get_enumeration_name(
            ctx, attribute);

        std::optional<Enumeration>
            enumeration = enumeration_name ?
                              std::make_optional(
                                  ArrayExperimental::get_enumeration(
                                      ctx, array, attribute.name())) :
                              std::nullopt;

        return std::make_shared<SOMABinaryColumn>(attribute, enumeration);

    } else if (soma_schema.contains(TILEDB_SOMA_SCHEMA_COL_DIM_KEY)) {
        std::vector<std::string>
            dimension_names = soma_schema[TILEDB_SOMA_SCHEMA_COL_DIM_KEY]
                                  .template get<std::vector<std::string>>();

        if (dimension_names.size() != 1) {
            throw TileDBSOMAError(fmt::format(
                "[SOMABinaryColumn][deserialize] Invalid number of dimensions: "
                "expected 1, got {}",
                dimension_names.size()));
        }

        auto dimension = array.schema().domain().dimension(dimension_names[0]);

        return std::make_shared<SOMABinaryColumn>(dimension);
    } else {
        throw TileDBSOMAError(
            "[SOMABinaryColumn][deserialize] Missing required field "
            "'tiledb_attributes' or 'tiledb_dimensions'");
    }
}

std::shared_ptr<SOMABinaryColumn> SOMABinaryColumn::create(
    std::shared_ptr<Context> ctx,
    ArrowSchema* schema,
    ArrowArray* array,
    const std::string& soma_type,
    bool is_index,
    std::string_view type_metadata,
    PlatformConfig platform_config) {
    if (is_index) {
        auto dimension = ArrowAdapter::tiledb_dimension_from_arrow_schema(
            ctx,
            schema,
            array,
            soma_type,
            type_metadata,
            "",
            "",
            platform_config);

        return std::make_shared<SOMABinaryColumn>(dimension);

    } else {
        auto attribute = ArrowAdapter::tiledb_attribute_from_arrow_schema(
            ctx, schema, type_metadata, platform_config);

        return std::make_shared<SOMABinaryColumn>(
            SOMABinaryColumn(attribute.first, attribute.second));
    }
}

void SOMABinaryColumn::_set_dim_points(
    ManagedQuery& query, const std::any& points) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_set_dim_points] Column with name {} is not an "
            "index "
            "column",
            name()));
    }

    std::vector<std::string> casted_points;

    for (const auto& bytes :
         std::any_cast<std::span<const std::vector<std::byte>>>(points)) {
        casted_points.push_back(std::string(
            reinterpret_cast<const char*>(bytes.data()), bytes.size()));
    }

    query.select_points(name(), casted_points);
}

void SOMABinaryColumn::_set_dim_ranges(
    ManagedQuery& query, const std::any& ranges) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_set_dim_ranges] Column with name {} is not an "
            "index "
            "column",
            name()));
    }

    std::vector<std::pair<std::string, std::string>> casted_ranges;

    for (const auto& bytes : std::any_cast<std::vector<
             std::pair<std::vector<std::byte>, std::vector<std::byte>>>>(
             ranges)) {
        casted_ranges.push_back(std::make_pair(
            std::string(
                reinterpret_cast<const char*>(bytes.first.data()),
                bytes.first.size()),
            std::string(
                reinterpret_cast<const char*>(bytes.second.data()),
                bytes.second.size())));
    }

    query.select_ranges(name(), casted_ranges);
}

void SOMABinaryColumn::_set_current_domain_slot(
    NDRectangle& rectangle, std::span<const std::any> domain) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_set_current_domain_slot] Column with name {} "
            "is "
            "not "
            "an index column",
            name()));
    }

    if (domain.size() != 1) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_set_current_domain_slot] Invalid domain size. "
            "Expected 1, got {}",
            domain.size()));
    }

    auto dom = std::any_cast<std::array<std::vector<std::byte>, 2>>(domain[0]);
    if (dom[0].size() == 0 && dom[1].size() == 0) {
        rectangle.set_range(name(), "", "\x7f");
    } else {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_set_current_domain_slot] domain (\"{}\", "
            "\"{}\") cannot be set for "
            "binary index columns: please use "
            "(\"\", \"\")",
            std::string(
                reinterpret_cast<const char*>(dom[0].data()), dom[0].size()),
            std::string(
                reinterpret_cast<const char*>(dom[1].data()), dom[1].size())));
    }
}

std::pair<bool, std::string> SOMABinaryColumn::_can_set_current_domain_slot(
    std::optional<NDRectangle>&, std::span<const std::any> new_domain) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_set_current_domain_slot] Column with name {} "
            "is "
            "not "
            "an index column",
            name()));
    }

    if (new_domain.size() != 1) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_can_set_current_domain_slot] Expected domain "
            "size for '{}' is 1, found {}",
            name(),
            new_domain.size()));
    }

    auto dom = std::any_cast<std::array<std::vector<std::byte>, 2>>(
        new_domain[0]);
    if (dom[0].size() != 0 || dom[1].size() != 0) {
        return std::pair(
            false,
            "domain cannot be set for string index columns: please use "
            "(\"\", \"\")");
    }

    return std::pair(true, "");
};

std::any SOMABinaryColumn::_core_domain_slot() const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_core_domain_slot] Column with name {} is not "
            "an "
            "index column",
            name()));
    }

    return std::make_any<
        std::pair<std::vector<std::byte>, std::vector<std::byte>>>(
        std::vector<std::byte>(), std::vector<std::byte>());
}

std::any SOMABinaryColumn::_non_empty_domain_slot(Array& array) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_non_empty_domain_slot] Column with name {} is "
            "not "
            "an "
            "index column",
            name()));
    }

    auto string_non_empty_domain = array.non_empty_domain_var(name());

    auto min_data_ptr = reinterpret_cast<std::byte*>(
        string_non_empty_domain.first.data());
    auto max_data_ptr = reinterpret_cast<std::byte*>(
        string_non_empty_domain.second.data());

    std::vector<std::byte> min(
        min_data_ptr, min_data_ptr + string_non_empty_domain.first.size());
    std::vector<std::byte> max(
        max_data_ptr, max_data_ptr + string_non_empty_domain.second.size());

    return std::make_any<
        std::pair<std::vector<std::byte>, std::vector<std::byte>>>(min, max);
}

std::any SOMABinaryColumn::_non_empty_domain_slot_opt(
    const SOMAContext& ctx, Array& array) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_non_empty_domain_slot] Column with name {} is "
            "not "
            "an "
            "index column",
            name()));
    }

    int32_t is_empty;
    void* var_start;
    void* var_end;
    uint64_t size_start, size_end;

    ctx.tiledb_ctx()->handle_error(
        tiledb_array_get_non_empty_domain_var_size_from_name(
            ctx.tiledb_ctx()->ptr().get(),
            array.ptr().get(),
            name().c_str(),
            &size_start,
            &size_end,
            &is_empty));

    if (is_empty) {
        return std::make_any<std::optional<
            std::pair<std::vector<std::byte>, std::vector<std::byte>>>>(
            std::nullopt);
    }

    var_start = malloc(size_start);
    var_end = malloc(size_end);

    ctx.tiledb_ctx()->handle_error(
        tiledb_array_get_non_empty_domain_var_from_name(
            ctx.tiledb_ctx()->ptr().get(),
            array.ptr().get(),
            name().c_str(),
            var_start,
            var_end,
            &is_empty));

    auto min_data_ptr = reinterpret_cast<std::byte*>(var_start);
    auto max_data_ptr = reinterpret_cast<std::byte*>(var_end);

    auto ned = std::make_pair(
        std::vector(min_data_ptr, min_data_ptr + size_start),
        std::vector(max_data_ptr, max_data_ptr + size_end));
    free(var_start);
    free(var_end);

    return std::make_any<std::optional<
        std::pair<std::vector<std::byte>, std::vector<std::byte>>>>(ned);
}

std::any SOMABinaryColumn::_core_current_domain_slot(
    const SOMAContext& ctx, Array& array) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_core_current_domain_slot] Column with name {} "
            "is "
            "not "
            "an index column",
            name()));
    }

    CurrentDomain
        current_domain = tiledb::ArraySchemaExperimental::current_domain(
            *ctx.tiledb_ctx(), array.schema());
    NDRectangle ndrect = current_domain.ndrectangle();

    return _core_current_domain_slot(ndrect);
}

std::any SOMABinaryColumn::_core_current_domain_slot(
    NDRectangle& ndrect) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][_core_current_domain_slot] Column with name {} "
            "is "
            "not "
            "an index column",
            name()));
    }

    std::array<std::string, 2> domain = ndrect.range<std::string>(name());

    if (domain[0] == "" && (domain[1] == "\x7f" || domain[1] == "\xff")) {
        return std::make_any<
            std::pair<std::vector<std::byte>, std::vector<std::byte>>>(
            std::vector<std::byte>(), std::vector<std::byte>());
    } else {
        throw TileDBSOMAError(fmt::format(
            "[SOMAColumn][core_current_domain_slot] unexpected current "
            "domain returnd ({}, {})",
            domain[0],
            domain[1]));
    }
}

std::pair<ArrowArray*, ArrowSchema*> SOMABinaryColumn::arrow_domain_slot(
    const SOMAContext& ctx, Array& array, enum Domainish kind) const {
    if (!isIndexColumn()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMABinaryColumn][arrow_domain_slot] Column with name {} is not "
            "an "
            "index column",
            name()));
    }

    ArrowArray* arrow_array = ArrowAdapter::make_arrow_array_child_binary(
        domain_slot<std::vector<std::byte>>(ctx, array, kind));

    return std::make_pair(arrow_array, arrow_schema_slot(ctx, array));
}

ArrowSchema* SOMABinaryColumn::arrow_schema_slot(
    const SOMAContext& ctx, Array& array) const {
    if (isIndexColumn()) {
        ArrowSchema* schema = ArrowAdapter::arrow_schema_from_tiledb_dimension(
            std::get<Dimension>(container));
        free((void*)schema->format);
        schema->format = strdup("Z");

        return schema;
    } else {
        return ArrowAdapter::arrow_schema_from_tiledb_attribute(
            std::get<Attribute>(container), *ctx.tiledb_ctx(), array);
    }
}

void SOMABinaryColumn::serialize(nlohmann::json& columns_schema) const {
    nlohmann::json column;

    column[TILEDB_SOMA_SCHEMA_COL_TYPE_KEY] = static_cast<uint32_t>(
        soma_column_datatype_t::SOMA_COLUMN_BINARY);

    if (isIndexColumn()) {
        column[TILEDB_SOMA_SCHEMA_COL_DIM_KEY] = {name()};
    } else {
        column[TILEDB_SOMA_SCHEMA_COL_ATTR_KEY] = {name()};
    }

    columns_schema.push_back(column);
}
}  // namespace tiledbsoma
