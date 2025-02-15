/**
 * @file   soma_column.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAColumn class.
 */

#include "soma_column.h"

#include "soma_attribute.h"
#include "soma_dimension.h"
#include "soma_geometry_column.h"

namespace tiledbsoma {

std::map<uint32_t, SOMAColumn::Factory> SOMAColumn::deserialiser_map = {
    {soma_column_datatype_t::SOMA_COLUMN_ATTRIBUTE, SOMAAttribute::deserialize},
    {soma_column_datatype_t::SOMA_COLUMN_DIMENSION, SOMADimension::deserialize},
    {soma_column_datatype_t::SOMA_COLUMN_GEOMETRY,
     SOMAGeometryColumn::deserialize}};

std::vector<std::shared_ptr<SOMAColumn>> SOMAColumn::deserialize(
    const Context& ctx,
    const Array& array,
    const std::map<std::string, tiledbsoma::MetadataValue>& metadata) {
    std::vector<std::shared_ptr<SOMAColumn>> columns;

    nlohmann::json soma_schema_columns = nlohmann::json::array();

    if (metadata.contains(TILEDB_SOMA_SCHEMA_KEY)) {
        auto soma_schema_extension_raw = metadata.at(TILEDB_SOMA_SCHEMA_KEY);
        auto data = static_cast<const char*>(
            std::get<2>(soma_schema_extension_raw));
        auto soma_schema_extension = data != nullptr ?
                                         nlohmann::json::parse(std::string(
                                             data,
                                             std::get<1>(
                                                 soma_schema_extension_raw))) :
                                         nlohmann::json::object();

        if (!soma_schema_extension.contains(TILEDB_SOMA_SCHEMA_COL_KEY)) {
            throw TileDBSOMAError(std::format(
                "[SOMAArray][fill_columns] Missing '{}' key from '{}'",
                TILEDB_SOMA_SCHEMA_COL_KEY,
                TILEDB_SOMA_SCHEMA_KEY));
        }

        soma_schema_columns = soma_schema_extension.value(
            TILEDB_SOMA_SCHEMA_COL_KEY, nlohmann::json::array());
    }

    if (!soma_schema_columns.empty()) {
        for (auto& column : soma_schema_columns) {
            auto type = column[TILEDB_SOMA_SCHEMA_COL_TYPE_KEY]
                            .template get<uint32_t>();

            auto col = deserialiser_map[type](column, ctx, array, metadata);

            if (col) {
                // Deserialized column can be null in case the array is modified
                // and the column no longer exists.
                columns.push_back(col);
            }
        }

        // Check for any newly added attributes
        std::unordered_set<std::string> used_attribute_names;

        std::for_each(
            columns.cbegin(),
            columns.cend(),
            [&used_attribute_names](const std::shared_ptr<SOMAColumn>& col) {
                if (col->tiledb_attributes().has_value()) {
                    auto attributes = col->tiledb_attributes().value();
                    for (const auto& attribute : attributes) {
                        used_attribute_names.insert(attribute.name());
                    }
                }
            });

        for (size_t i = 0; i < array.schema().attribute_num(); ++i) {
            auto attribute = array.schema().attribute(i);

            // Attribute is already used by another attribute so we skip
            if (used_attribute_names.contains(attribute.name())) {
                continue;
            }

            auto enumeration_name = AttributeExperimental::get_enumeration_name(
                ctx, attribute);
            auto enumeration = enumeration_name.has_value() ?
                                   std::make_optional(
                                       ArrayExperimental::get_enumeration(
                                           ctx,
                                           array,
                                           enumeration_name.value())) :
                                   std::nullopt;

            columns.push_back(
                std::make_shared<SOMAAttribute>(attribute, enumeration));
        }
    } else {
        // All arrays before the introduction of SOMAColumn do not have
        // composite columns, thus the metadata are trivially constructible
        for (auto& dimension : array.schema().domain().dimensions()) {
            columns.push_back(std::make_shared<SOMADimension>(dimension));
        }

        for (size_t i = 0; i < array.schema().attribute_num(); ++i) {
            auto attribute = array.schema().attribute(i);
            auto enumeration_name = AttributeExperimental::get_enumeration_name(
                ctx, attribute);
            auto enumeration = enumeration_name.has_value() ?
                                   std::make_optional(
                                       ArrayExperimental::get_enumeration(
                                           ctx,
                                           array,
                                           enumeration_name.value())) :
                                   std::nullopt;

            columns.push_back(
                std::make_shared<SOMAAttribute>(attribute, enumeration));
        }
    }

    return columns;
}

template <>
std::pair<std::string, std::string> SOMAColumn::core_domain_slot<std::string>()
    const {
    return std::pair<std::string, std::string>("", "");
}

template <>
std::pair<std::string, std::string>
SOMAColumn::core_current_domain_slot<std::string>(
    const SOMAContext& ctx, Array& array) const {
    // Here is an intersection of a few oddities:
    //
    // * Core domain for string dims must be a nullptr pair; it cannot
    //   be anything else.
    // * TileDB-Py shows this by using an empty-string pair, which we
    //   imitate.
    // * Core current domain for string dims must _not_ be a nullptr
    //   pair.
    // * In TileDB-SOMA, unless the user specifies otherwise, we use ""
    //   for min and "\x7f" for max. (We could use "\x7f" but that causes
    //   display problems in Python.)
    //
    // To work with all these factors, if the current domain is the default ""
    // to "\x7f", return an empty-string pair just as we do for domain. (There
    // was some pre-1.15 software using "\xff" and it's super-cheap to check for
    // that as well.)
    try {
        std::pair<std::string, std::string>
            current_domain = std::any_cast<std::pair<std::string, std::string>>(
                _core_current_domain_slot(ctx, array));

        if (current_domain.first == "" && (current_domain.second == "\x7f" ||
                                           current_domain.second == "\xff")) {
            return std::pair<std::string, std::string>("", "");
        } else {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][core_current_domain_slot] unexpected current "
                "domain returnd ({}, {})",
                current_domain.first,
                current_domain.second));
        }
    } catch (const std::exception& e) {
        throw TileDBSOMAError(e.what());
    }
}

template <>
std::pair<std::string, std::string>
SOMAColumn::core_current_domain_slot<std::string>(NDRectangle& ndrect) const {
    try {
        std::pair<std::string, std::string>
            current_domain = std::any_cast<std::pair<std::string, std::string>>(
                _core_current_domain_slot(ndrect));

        if (current_domain.first == "" && (current_domain.second == "\x7f" ||
                                           current_domain.second == "\xff")) {
            return std::pair<std::string, std::string>("", "");
        } else {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][core_current_domain_slot] unexpected current "
                "domain returnd ({}, {})",
                current_domain.first,
                current_domain.second));
        }
    } catch (const std::exception& e) {
        throw TileDBSOMAError(e.what());
    }
}
}  // namespace tiledbsoma
