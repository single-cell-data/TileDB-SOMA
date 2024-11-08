#ifndef SOMA_ATTRIBUTE
#define SOMA_ATTRIBUTE

#include <algorithm>
#include <vector>

#include <tiledb/tiledb>
#include "soma_column.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAAttribute : public virtual SOMAColumn {
   public:
    SOMAAttribute(const Context& ctx, Array& array, Attribute attribute)
        : SOMAColumn(ctx, array)
        , attribute(attribute) {
    }

    inline std::string domain_to_str() const {
        return "";
    }

    virtual inline std::string name() const {
        return attribute.name();
    }

    virtual inline bool isIndex() const {
        return false;
    }

   private:
    virtual void _set_dim_ranges(
        ManagedQuery& query, const std::any& ranges) const;

    virtual std::any _core_domain_slot() const;

    virtual std::any _non_empty_domain_slot() const;

    virtual std::any _core_current_domain_slot() const;

    Attribute attribute;
};
}  // namespace tiledbsoma

#endif