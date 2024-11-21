#ifndef SOMA_ATTRIBUTE_H
#define SOMA_ATTRIBUTE_H

#include <algorithm>
#include <vector>

#include <tiledb/tiledb>
#include "soma_column.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAAttribute : public virtual SOMAColumn {
   public:
    /**
     * Create a ``SOMAAttribute`` shared pointer from an arrow schema
     */
    static std::shared_ptr<SOMAAttribute> create(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        std::string_view type_metadata,
        PlatformConfig platform_config);

    SOMAAttribute(
        Attribute attribute,
        std::optional<Enumeration> enumeration = std::nullopt)
        : attribute(attribute)
        , enumeration(enumeration) {
    }

    virtual inline std::string name() const {
        return attribute.name();
    }

    virtual inline bool isIndexColumn() const {
        return false;
    }

    virtual inline void select_columns(
        const std::unique_ptr<ManagedQuery>& query,
        bool if_not_empty = false) const override {
        query->select_columns(std::vector({attribute.name()}), if_not_empty);
    };

    inline soma_column_datatype_t type() const {
        return soma_column_datatype_t::SOMA_COLUMN_ATTRIBUTE;
    }

    inline std::optional<tiledb_datatype_t> domain_type() const {
        return std::nullopt;
    }

    inline std::optional<tiledb_datatype_t> data_type() const {
        return attribute.type();
    }

    inline std::optional<std::vector<Dimension>> tiledb_dimensions() {
        return std::nullopt;
    }

    inline std::optional<std::vector<Attribute>> tiledb_attributes() {
        return std::vector({attribute});
    }

    inline virtual std::optional<std::vector<Enumeration>>
    tiledb_enumerations() {
        if (!enumeration.has_value()) {
            return std::nullopt;
        }

        return std::vector({enumeration.value()});
    }

    virtual ArrowArray* arrow_domain_slot(
        const SOMAContext& ctx,
        Array& array,
        enum Domainish kind) const override;

    virtual ArrowSchema* arrow_schema_slot(
        const SOMAContext& ctx, Array& array) override;

   private:
    virtual void _set_dim_points(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::any& points) const;

    virtual void _set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::any& ranges) const;

    virtual void _set_current_domain_slot(
        NDRectangle& rectangle, const std::vector<const void*>& domain) const;

    virtual std::any _core_domain_slot() const;

    virtual std::any _non_empty_domain_slot(Array& array) const;

    virtual std::any _core_current_domain_slot(
        const SOMAContext& ctx, Array& array) const;

    Attribute attribute;
    std::optional<Enumeration> enumeration;
};
}  // namespace tiledbsoma

#endif