/**
 * @file   soma_attribute.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAAttribute class. SOMAAttribute extends SOMAColumn
 * and wraps a TileDB Attribute and optionally an enumeration associated with
 * the attribute. The purpose of this class is to provide a common interface
 * identical to TileDB dimensions and composite columns.
 */

#ifndef SOMA_ATTRIBUTE_H
#define SOMA_ATTRIBUTE_H

#include <algorithm>
#include <vector>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include "../tiledb_adapter/platform_config.h"
#include "soma_column.h"

namespace tiledbsoma {

class SOMAAttribute : public SOMAColumn {
   public:
    //===================================================================
    //= public static
    //===================================================================

    static std::shared_ptr<SOMAColumn> deserialize(
        const nlohmann::json& soma_schema,
        const tiledb::Context& ctx,
        const tiledb::Array& array,
        const std::map<std::string, tiledbsoma::MetadataValue>& metadata);

    /**
     * Create a ``SOMAAttribute`` shared pointer from an Arrow schema
     */
    static std::shared_ptr<SOMAAttribute> create(
        std::shared_ptr<tiledb::Context> ctx,
        ArrowSchema* schema,
        std::string_view type_metadata,
        PlatformConfig platform_config);

    SOMAAttribute(tiledb::Attribute attribute, std::optional<tiledb::Enumeration> enumeration = std::nullopt)
        : attribute(attribute)
        , enumeration(enumeration) {
    }

    inline std::string name() const override {
        return attribute.name();
    }

    inline bool isIndexColumn() const override {
        return false;
    }

    inline void select_columns(common::ManagedQuery& query, bool if_not_empty = false) const override {
        query.select_columns(std::vector({attribute.name()}), if_not_empty);
    };

    inline soma_column_datatype_t type() const override {
        return soma_column_datatype_t::SOMA_COLUMN_ATTRIBUTE;
    }

    inline std::optional<tiledb_datatype_t> domain_type() const override {
        return std::nullopt;
    }

    inline std::optional<tiledb_datatype_t> data_type() const override {
        return attribute.type();
    }

    inline std::optional<std::vector<tiledb::Dimension>> tiledb_dimensions() override {
        return std::nullopt;
    }

    inline std::optional<std::vector<tiledb::Attribute>> tiledb_attributes() override {
        return std::vector({attribute});
    }

    inline std::optional<std::vector<tiledb::Enumeration>> tiledb_enumerations() override {
        if (!enumeration.has_value()) {
            return std::nullopt;
        }

        return std::vector({enumeration.value()});
    }

    std::pair<ArrowArray*, ArrowSchema*> arrow_domain_slot(
        const SOMAContext& ctx,
        tiledb::Array& array,
        enum Domainish kind,
        bool downcast_dict_of_large_var = false) const override;

    ArrowSchema* arrow_schema_slot(
        const SOMAContext& ctx, tiledb::Array& array, bool downcast_dict_of_large_var = false) const override;

    void serialize(nlohmann::json&) const override;

   private:
    void _set_dim_points(common::ManagedQuery& query, const std::any& points) const override;

    void _set_dim_ranges(common::ManagedQuery& query, const std::any& ranges) const override;

    void _set_current_domain_slot(tiledb::NDRectangle& rectangle, std::span<const std::any> domain) const override;

    std::pair<bool, std::string> _can_set_current_domain_slot(
        std::optional<tiledb::NDRectangle>& rectangle, std::span<const std::any> new_domain) const override;

    std::any _core_domain_slot() const override;

    std::any _non_empty_domain_slot(tiledb::Array& array) const override;

    std::any _non_empty_domain_slot_opt(const SOMAContext& ctx, tiledb::Array& array) const override;

    std::any _core_current_domain_slot(const SOMAContext& ctx, tiledb::Array& array) const override;

    std::any _core_current_domain_slot(tiledb::NDRectangle& ndrect) const override;

    tiledb::Attribute attribute;
    std::optional<tiledb::Enumeration> enumeration;
};
}  // namespace tiledbsoma

#endif
