/**
 * @file   soma_geometry_column.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAGeometryColumn class. SOMAGeometryColumn wraps a
 * TileDB Attribute storing the WKB (Well-Known Binary) encoded geometry and
 * adds a collection of internal TileDB dimension to provide spatial indexing.
 * It implements function to perform queries as well as core domain and current
 * domain operations. The purpose of this class is to provide a common interface
 * identical to TileDB dimensions, attributes and other composite columns.
 *
 * The current indexing mechanish adapts the idea of Priority R-tree on top of a
 * TileDB Array using a set of TIleDB Dimensions to store the MBR corners of
 * each geometry.
 */

#ifndef SOMA_GEOMETRY_COLUMN_H
#define SOMA_GEOMETRY_COLUMN_H

#include <algorithm>
#include <vector>

#include <tiledb/tiledb>
#include "soma_column.h"
#include "soma_coordinates.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMAGeometryColumn : public SOMAColumn {
   public:
    //===================================================================
    //= public static
    //===================================================================
    static std::shared_ptr<SOMAColumn> deserialize(
        const nlohmann::json& soma_schema,
        const Context& ctx,
        const Array& array,
        const std::map<std::string, tiledbsoma::MetadataValue>&);

    static std::shared_ptr<SOMAGeometryColumn> create(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        ArrowSchema* spatial_schema,
        ArrowArray* spatial_array,
        const SOMACoordinateSpace& coordinate_space,
        const std::string& soma_type,
        std::string_view type_metadata,
        PlatformConfig platform_config);

    SOMAGeometryColumn(
        std::vector<Dimension> dimensions,
        Attribute attribute,
        SOMACoordinateSpace coordinate_space)
        : dimensions(dimensions)
        , attribute(attribute)
        , coordinate_space(coordinate_space){};

    inline std::string name() const override {
        return SOMA_GEOMETRY_COLUMN_NAME;
    }

    inline bool isIndexColumn() const override {
        return true;
    }

    inline void select_columns(
        ManagedQuery& query, bool if_not_empty = false) const override {
        query.select_columns(std::vector({attribute.name()}), if_not_empty);
    };

    inline soma_column_datatype_t type() const override {
        return soma_column_datatype_t::SOMA_COLUMN_GEOMETRY;
    }

    inline std::optional<tiledb_datatype_t> domain_type() const override {
        return dimensions.front().type();
    }

    inline std::optional<tiledb_datatype_t> data_type() const override {
        return attribute.type();
    }

    inline std::optional<std::vector<Dimension>> tiledb_dimensions() override {
        return dimensions;
    }

    inline std::optional<std::vector<Attribute>> tiledb_attributes() override {
        return std::vector({attribute});
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

    inline SOMACoordinateSpace get_coordinate_space() const {
        return coordinate_space;
    }

   protected:
    void _set_dim_points(
        ManagedQuery& query, const std::any& points) const override;

    void _set_dim_ranges(
        ManagedQuery& query, const std::any& ranges) const override;

    void _set_current_domain_slot(
        NDRectangle& rectangle,
        std::span<const std::any> new_current_domain) const override;

    std::pair<bool, std::string> _can_set_current_domain_slot(
        std::optional<NDRectangle>& rectangle,
        std::span<const std::any> new_current_domain) const override;

    std::any _core_domain_slot() const override;

    std::any _non_empty_domain_slot(Array& array) const override;

    std::any _non_empty_domain_slot_opt(
        const SOMAContext& ctx, Array& array) const override;

    std::any _core_current_domain_slot(
        const SOMAContext& ctx, Array& array) const override;

    std::any _core_current_domain_slot(NDRectangle& ndrect) const override;

   private:
    /**
     * The current implementation of SOMAGeometryColumn uses a pair of TileDB
     * dimensions to store the min and max point of the bounding box per
     * dimension. E.g. a 2D geometry will have 4 TileDB dimensions (2 *
     * num_spatial_axes) to provide spatial indexing.
     */
    const size_t TDB_DIM_PER_SPATIAL_AXIS = 2;
    std::vector<Dimension> dimensions;
    Attribute attribute;
    SOMACoordinateSpace coordinate_space;

    /**
     * Compute the usable domain limits. If the array has a current domain then
     * it is used to compute the limits, otherwise the core domain is used.
     */
    std::vector<std::pair<double_t, double_t>> _limits(
        const Context& ctx, const ArraySchema& schema) const;

    std::vector<std::pair<double_t, double_t>> _transform_ranges(
        const std::vector<
            std::pair<std::vector<double_t>, std::vector<double_t>>>& ranges)
        const;

    std::vector<std::pair<double_t, double_t>> _transform_points(
        const std::span<const std::vector<double_t>>& points) const;
};

}  // namespace tiledbsoma
#endif
