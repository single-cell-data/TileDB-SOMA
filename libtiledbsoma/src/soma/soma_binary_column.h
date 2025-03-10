/**
 * @file   soma_binary_column.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMABinaryColumn class. SOMABinaryColumn acts as a
 * wrapper to a TileDB Dimension holding binary data encoded as string and
 * implements function to perform queries as well as core domain and current
 * domain operations. It provides a common interface identical to TileDB
 * attributes, dimensions and composite columns.
 */

#ifndef SOMA_BINARY_COLUMN_H
#define SOMA_BINARY_COLUMN_H

#include <algorithm>
#include <variant>
#include <vector>

#include <tiledb/tiledb>
#include "soma_column.h"

namespace tiledbsoma {

using namespace tiledb;

class SOMABinaryColumn : public SOMAColumn {
   public:
    //===================================================================
    //= public static
    //===================================================================

    static std::shared_ptr<SOMAColumn> deserialize(
        const nlohmann::json& soma_schema,
        const Context& ctx,
        const Array& array,
        const std::map<std::string, tiledbsoma::MetadataValue>& metadata);

    static std::shared_ptr<SOMABinaryColumn> create(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        ArrowArray* array,
        const std::string& soma_type,
        bool is_index,
        std::string_view type_metadata,
        PlatformConfig platform_config);

    SOMABinaryColumn(Dimension dimension)
        : container(dimension) {
    }

    SOMABinaryColumn(
        Attribute attribute, std::optional<Enumeration> enumeration)
        : container(attribute)
        , enumeration(enumeration) {
    }

    inline std::string name() const override {
        return std::visit([](auto&& x) { return x.name(); }, container);
    }

    inline bool isIndexColumn() const override {
        return std::holds_alternative<Dimension>(container);
    }

    inline void select_columns(
        ManagedQuery& query, bool if_not_empty = false) const override {
        query.select_columns(std::vector({name()}), if_not_empty);
    };

    inline soma_column_datatype_t type() const override {
        return soma_column_datatype_t::SOMA_COLUMN_BINARY;
    }

    inline std::optional<tiledb_datatype_t> domain_type() const override {
        if (std::holds_alternative<Dimension>(container)) {
            return TILEDB_BLOB;
        }

        return std::nullopt;
    }

    inline std::optional<tiledb_datatype_t> data_type() const override {
        if (std::holds_alternative<Attribute>(container)) {
            return TILEDB_BLOB;
        }

        return std::nullopt;
    }

    inline std::optional<std::vector<Dimension>> tiledb_dimensions() override {
        if (std::holds_alternative<Dimension>(container)) {
            return std::vector({std::get<Dimension>(container)});
        }

        return std::nullopt;
    }

    inline std::optional<std::vector<Attribute>> tiledb_attributes() override {
        if (std::holds_alternative<Attribute>(container)) {
            return std::vector({std::get<Attribute>(container)});
        }

        return std::nullopt;
    }

    inline std::optional<std::vector<Enumeration>> tiledb_enumerations()
        override {
        if (!enumeration.has_value()) {
            return std::nullopt;
        }

        return std::vector({enumeration.value()});
    }

    std::pair<ArrowArray*, ArrowSchema*> arrow_domain_slot(
        const SOMAContext& ctx,
        Array& array,
        enum Domainish kind) const override;

    ArrowSchema* arrow_schema_slot(
        const SOMAContext& ctx, Array& array) const override;

    void serialize(nlohmann::json&) const override;

   protected:
    void _set_dim_points(
        ManagedQuery& query, const std::any& ranges) const override;

    void _set_dim_ranges(
        ManagedQuery& query, const std::any& ranges) const override;

    void _set_current_domain_slot(
        NDRectangle& rectangle,
        std::span<const std::any> domain) const override;

    std::pair<bool, std::string> _can_set_current_domain_slot(
        std::optional<NDRectangle>& rectangle,
        std::span<const std::any> new_domain) const override;

    std::any _core_domain_slot() const override;

    std::any _non_empty_domain_slot(Array& array) const override;

    std::any _non_empty_domain_slot_opt(
        const SOMAContext& ctx, Array& array) const override;

    std::any _core_current_domain_slot(
        const SOMAContext& ctx, Array& array) const override;

    std::any _core_current_domain_slot(NDRectangle& ndrect) const override;

   private:
    std::variant<Attribute, Dimension> container;
    std::optional<Enumeration> enumeration;
};
}  // namespace tiledbsoma

#endif