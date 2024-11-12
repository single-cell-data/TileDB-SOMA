#ifndef SOMA_DIMENSION
#define SOMA_DIMENSION

#include <algorithm>
#include <vector>

#include <tiledb/tiledb>
#include "soma_column.h"

namespace tiledbsoma {

using namespace tiledb;

class SOMADimension : public virtual SOMAColumn {
   public:
    static std::shared_ptr<SOMADimension> create(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        ArrowArray* array,
        const std::string& soma_type,
        std::string_view type_metadata,
        PlatformConfig platform_config);

    SOMADimension(const Context& ctx, Dimension dimension)
        : SOMAColumn(ctx)
        , dimension(dimension) {
    }

    inline std::string name() const {
        return dimension.name();
    }

    inline bool isIndex() const {
        return true;
    }

    inline std::string domain_to_str() const {
        return "";
    }

    inline std::optional<std::vector<Dimension>> tiledb_dimensions() {
        return std::vector({dimension});
    }

    inline std::optional<std::vector<Attribute>> tiledb_attributes() {
        return std::nullopt;
    }

    inline virtual std::optional<std::vector<Enumeration>>
    tiledb_enumerations() {
        return std::nullopt;
    }

    inline virtual bool has_current_domain() {
        return _has_current_domain;
    }

   protected:
    virtual void _set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const std::any& ranges) const;

    virtual void _set_current_domain_slot(
        NDRectangle& rectangle, const std::vector<const void*>& domain) const;

    virtual std::any _core_domain_slot() const;

    virtual std::any _non_empty_domain_slot(Array& array) const;

    virtual std::any _core_current_domain_slot(Array& array) const;

   private:
    Dimension dimension;

    bool _has_current_domain;
};
}  // namespace tiledbsoma

#endif