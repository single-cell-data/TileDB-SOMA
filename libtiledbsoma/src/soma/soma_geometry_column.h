#ifndef SOMA_GEOMETRY_COLUMN_H
#define SOMA_GEOMETRY_COLUMN_H

#include <algorithm>
#include <vector>

#include <tiledb/tiledb>
#include "soma_column.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMAGeometryColumn : public virtual SOMAColumn {
   public:
    static std::shared_ptr<SOMAGeometryColumn> create(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        ArrowArray* array,
        ArrowSchema* spatial_schema,
        ArrowArray* spatial_array,
        const std::string& soma_type,
        std::string_view type_metadata,
        PlatformConfig platform_config);

    SOMAGeometryColumn(std::vector<Dimension> dimensions, Attribute attribute)
        : dimensions(dimensions)
        , attribute(attribute){};

    virtual inline std::string name() const override {
        return SOMA_GEOMETRY_COLUMN_NAME;
    }

    virtual inline bool isIndexColumn() const override {
        return true;
    }

    virtual inline void select_columns(
        const std::unique_ptr<ManagedQuery>& query,
        bool if_not_empty = false) const override {
        query->select_columns(std::vector({attribute.name()}), if_not_empty);
    };

    virtual inline soma_column_datatype_t type() const override {
        return soma_column_datatype_t::SOMA_COLUMN_GEOMETRY;
    }

    virtual inline std::optional<tiledb_datatype_t> domain_type()
        const override {
        return dimensions.front().type();
    }

    virtual inline std::optional<tiledb_datatype_t> data_type() const override {
        return attribute.type();
    }

    virtual inline std::optional<std::vector<Dimension>> tiledb_dimensions()
        override {
        return dimensions;
    }

    virtual inline std::optional<std::vector<Attribute>> tiledb_attributes()
        override {
        return std::vector({attribute});
    }

    virtual inline std::optional<std::vector<Enumeration>> tiledb_enumerations()
        override {
        return std::nullopt;
    }

    virtual ArrowArray* arrow_domain_slot(
        const SOMAContext& ctx,
        Array& array,
        enum Domainish kind) const override;

    virtual ArrowSchema* arrow_schema_slot(
        const SOMAContext& ctx, Array& array) override;

   protected:
    virtual void _set_dim_points(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::any& points) const override;

    virtual void _set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::any& ranges) const override;

    virtual void _set_current_domain_slot(
        NDRectangle& rectangle,
        std::span<const std::any> domain) const override;

    virtual std::any _core_domain_slot() const override;

    virtual std::any _non_empty_domain_slot(Array& array) const override;

    virtual std::any _core_current_domain_slot(
        const SOMAContext& ctx, Array& array) const override;

   private:
    std::vector<Dimension> dimensions;
    Attribute attribute;

    std::vector<std::pair<double_t, double_t>> _limits(
        const SOMAContext& ctx, const ArraySchema& schema) const;

    std::vector<std::pair<double_t, double_t>> _transform_ranges(
        const std::vector<
            std::pair<std::vector<double_t>, std::vector<double_t>>>& ranges)
        const;

    std::vector<std::pair<double_t, double_t>> _transform_points(
        const std::span<const std::vector<double_t>>& points) const;
};

}  // namespace tiledbsoma
#endif