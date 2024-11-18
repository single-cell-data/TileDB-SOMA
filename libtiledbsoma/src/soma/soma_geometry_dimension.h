#ifndef SOMA_GEOMETRY_COLUMN
#define SOMA_GEOMETRY_COLUMN

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

    SOMAGeometryColumn(
        const Context& ctx,
        std::vector<Dimension> dimensions,
        Attribute attribute)
        : SOMAColumn(ctx)
        , dimensions(dimensions)
        , attribute(attribute){};

    inline std::string name() const {
        return SOMA_GEOMETRY_COLUMN_NAME;
    }

    inline std::string domain_to_str() const {
        return "";
    }

    inline bool isIndexColumn() const {
        return true;
    }

    inline virtual void select_columns(
        const std::unique_ptr<ManagedQuery>& query,
        bool if_not_empty = false) const override {
        query->select_columns(std::vector({attribute.name()}), if_not_empty);
    };

    inline soma_column_datatype_t type() const {
        return soma_column_datatype_t::SOMA_COLUMN_GEOMETRY;
    }

    inline std::optional<tiledb_datatype_t> domain_type() const {
        return dimensions.front().type();
    }

    inline std::optional<tiledb_datatype_t> data_type() const {
        return attribute.type();
    }

    inline std::optional<std::vector<Dimension>> tiledb_dimensions() {
        return dimensions;
    }

    inline std::optional<std::vector<Attribute>> tiledb_attributes() {
        return std::vector({attribute});
    }

    inline virtual std::optional<std::vector<Enumeration>>
    tiledb_enumerations() {
        return std::nullopt;
    }

    inline virtual bool has_current_domain() {
        return _has_current_domain;
    }

    virtual ArrowArray* arrow_domain_slot(
        Array& array, enum Domainish kind) const override;

    virtual ArrowSchema* arrow_schema_slot(
        const Context& ctx, Array& array) override;

   protected:
    virtual void _set_dim_points(
        const std::unique_ptr<ManagedQuery>& query,
        const std::any& points) const;

    void _set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const std::any& ranges) const;

    virtual void _set_current_domain_slot(
        NDRectangle& rectangle, const std::vector<const void*>& domain) const;

    virtual std::any _core_domain_slot() const;

    virtual std::any _non_empty_domain_slot(Array& array) const;

    virtual std::any _core_current_domain_slot(Array& array) const;

   private:
    std::vector<Dimension> dimensions;
    Attribute attribute;
    bool _has_current_domain;

    std::vector<std::pair<double_t, double_t>> _limits(
        const ArraySchema& schema) const;

    std::vector<std::pair<double_t, double_t>> _transform_ranges(
        const std::vector<
            std::pair<std::vector<double_t>, std::vector<double_t>>>& ranges)
        const;

    std::vector<std::pair<double_t, double_t>> _transform_points(
        const std::vector<std::vector<double_t>>& points) const;
};

}  // namespace tiledbsoma
#endif