/**
 * @file   soma_attribute.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAAttribute class.
 */

#include "soma_attribute.h"

namespace tiledbsoma {
std::shared_ptr<SOMAColumn> SOMAAttribute::deserialize(
    const nlohmann::json& soma_schema,
    const Context& ctx,
    const Array& array,
    const std::map<std::string, tiledbsoma::MetadataValue>&) {
    if (!soma_schema.contains(TILEDB_SOMA_SCHEMA_COL_ATTR_KEY)) {
        throw TileDBSOMAError(
            "[SOMAAttribute][deserialize] Missing required field "
            "'tiledb_attributes'");
    }

    std::vector<std::string>
        attribute_names = soma_schema[TILEDB_SOMA_SCHEMA_COL_ATTR_KEY]
                              .template get<std::vector<std::string>>();

    if (attribute_names.size() != 1) {
        throw TileDBSOMAError(std::format(
            "[SOMAAttribute][deserialize] Invalid number of attributes. "
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
        enumeration = enumeration_name.has_value() ?
                          std::make_optional(ArrayExperimental::get_enumeration(
                              ctx, array, enumeration_name.value())) :
                          std::nullopt;

    return std::make_shared<SOMAAttribute>(attribute, enumeration);
}

std::shared_ptr<SOMAAttribute> SOMAAttribute::create(
    std::shared_ptr<Context> ctx,
    ArrowSchema* schema,
    std::string_view type_metadata,
    PlatformConfig platform_config) {
    auto attribute = ArrowAdapter::tiledb_attribute_from_arrow_schema(
        ctx, schema, type_metadata, platform_config);

    return std::make_shared<SOMAAttribute>(
        SOMAAttribute(attribute.first, attribute.second));
}

void SOMAAttribute::_set_dim_points(ManagedQuery&, const std::any&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_set_dim_points] Column with name {} is not an index "
        "column",
        name()));
}

void SOMAAttribute::_set_dim_ranges(ManagedQuery&, const std::any&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_set_dim_ranges] Column with name {} is not an index "
        "column",
        name()));
}

void SOMAAttribute::_set_current_domain_slot(
    NDRectangle&, std::span<const std::any>) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_set_current_domain_slot] Column with name {} is not "
        "an index column",
        name()));
}

std::pair<bool, std::string> SOMAAttribute::_can_set_current_domain_slot(
    std::optional<NDRectangle>&, std::span<const std::any>) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_set_current_domain_slot] Column with name {} is not "
        "an index column",
        name()));
};

std::any SOMAAttribute::_core_domain_slot() const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_core_domain_slot] Column with name {} is not an "
        "index column",
        name()));
}

std::any SOMAAttribute::_non_empty_domain_slot(Array&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_non_empty_domain_slot] Column with name {} is not an "
        "index column",
        name()));
}

std::any SOMAAttribute::_non_empty_domain_slot_opt(
    const SOMAContext&, Array&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_non_empty_domain_slot] Column with name {} is not an "
        "index column",
        name()));
}

std::any SOMAAttribute::_core_current_domain_slot(
    const SOMAContext&, Array&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_core_current_domain_slot] Column with name {} is not "
        "an index column",
        name()));
}

std::any SOMAAttribute::_core_current_domain_slot(NDRectangle&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_core_current_domain_slot] Column with name {} is not "
        "an index column",
        name()));
}

std::pair<ArrowArray*, ArrowSchema*> SOMAAttribute::arrow_domain_slot(
    const SOMAContext&, Array&, enum Domainish) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][arrow_domain_slot] Column with name {} is not an "
        "index column",
        name()));
}

ArrowSchema* SOMAAttribute::arrow_schema_slot(
    const SOMAContext& ctx, Array& array) const {
    return ArrowAdapter::arrow_schema_from_tiledb_attribute(
               attribute, *ctx.tiledb_ctx(), array)
        .release();
}

void SOMAAttribute::serialize(nlohmann::json& columns_schema) const {
    nlohmann::json column;

    column[TILEDB_SOMA_SCHEMA_COL_TYPE_KEY] = static_cast<uint32_t>(
        soma_column_datatype_t::SOMA_COLUMN_ATTRIBUTE);
    column[TILEDB_SOMA_SCHEMA_COL_ATTR_KEY] = {attribute.name()};

    columns_schema.push_back(column);
}
}  // namespace tiledbsoma
