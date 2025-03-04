/**
 * @file   soma_dimension.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMADimension class. SOMADimension acts as a wrapper
 * to a TileDB Dimension and implements function to perform queries as well as
 * core domain and current domain operations. It provides a common interface
 * identical to TileDB attributes and composite columns.
 */

#ifndef SOMA_DIMENSION_H
#define SOMA_DIMENSION_H

#include <algorithm>
#include <vector>

#include <tiledb/tiledb>
#include "soma_column.h"

namespace tiledbsoma {

using namespace tiledb;

class SOMADimension : public SOMAColumn {
   public:
    //===================================================================
    //= public static
    //===================================================================

    static std::shared_ptr<SOMAColumn> deserialize(
        const nlohmann::json& soma_schema,
        const Context& ctx,
        const Array& array,
        const std::map<std::string, tiledbsoma::MetadataValue>& metadata);

    static std::shared_ptr<SOMADimension> create(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        ArrowArray* array,
        const std::string& soma_type,
        std::string_view type_metadata,
        PlatformConfig platform_config);

    SOMADimension(Dimension dimension)
        : dimension(dimension) {
    }

    inline std::string name() const override {
        return dimension.name();
    }

    inline bool isIndexColumn() const override {
        return true;
    }

    inline void select_columns(
        ManagedQuery& query, bool if_not_empty = false) const override {
        query.select_columns(std::vector({dimension.name()}), if_not_empty);
    };

    inline soma_column_datatype_t type() const override {
        return soma_column_datatype_t::SOMA_COLUMN_DIMENSION;
    }

    inline std::optional<tiledb_datatype_t> domain_type() const override {
        return dimension.type();
    }

    inline std::optional<tiledb_datatype_t> data_type() const override {
        return std::nullopt;
    }

    inline std::optional<std::vector<Dimension>> tiledb_dimensions() override {
        return std::vector({dimension});
    }

    inline std::optional<std::vector<Attribute>> tiledb_attributes() override {
        return std::nullopt;
    }

    inline std::optional<std::vector<Enumeration>> tiledb_enumerations()
        override {
        return std::nullopt;
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
    Dimension dimension;
};
}  // namespace tiledbsoma

#endif
