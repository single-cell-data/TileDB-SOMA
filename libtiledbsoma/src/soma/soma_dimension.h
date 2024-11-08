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
    SOMADimension(const Context& ctx, Array& array, Dimension dimension)
        : SOMAColumn(ctx, array)
        , dimension(dimension) {
    }

    inline std::string domain_to_str() const {
        return "";
    }

   protected:
    virtual void _set_dim_ranges(
        ManagedQuery& query, const std::any& ranges) const;

    virtual std::any _core_domain_slot() const;

    virtual std::any _non_empty_domain_slot() const;

    virtual std::any _core_current_domain_slot() const;

   private:
    Dimension dimension;
};
}  // namespace tiledbsoma

#endif