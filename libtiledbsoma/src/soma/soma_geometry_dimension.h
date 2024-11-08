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
    static SOMAGeometryColumn create(const Context& ctx, Array& array);

    SOMAGeometryColumn(
        const Context& ctx,
        Array& array,
        std::vector<Dimension> dimensions,
        Attribute attribute)
        : SOMAColumn(ctx, array)
        , dimensions(dimensions)
        , attribute(attribute){};

    inline std::string name() const {
        return SOMA_GEOMETRY_COLUMN_NAME;
    }

    inline std::string domain_to_str() const {
        return "";
    }

    inline bool isIndex() const {
        return true;
    }

   protected:
    void _set_dim_ranges(ManagedQuery& query, const std::any& ranges) const;

    virtual std::any _core_domain_slot() const;

    virtual std::any _non_empty_domain_slot() const;

    virtual std::any _core_current_domain_slot() const;

   private:
    std::vector<Dimension> dimensions;
    Attribute attribute;

    std::vector<std::pair<double_t, double_t>> _limits() const;

    std::vector<std::pair<double_t, double_t>> _transform_ranges(
        const std::vector<
            std::pair<std::vector<double_t>, std::vector<double_t>>>& ranges)
        const;
};

}  // namespace tiledbsoma
#endif