#include "soma_attribute.h"

namespace tiledbsoma {
void SOMAAttribute::_set_dim_ranges(
    ManagedQuery& query, const std::any& ranges) const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}

std::any SOMAAttribute::_core_domain_slot() const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}

std::any SOMAAttribute::_non_empty_domain_slot() const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}

std::any SOMAAttribute::_core_current_domain_slot() const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}
}  // namespace tiledbsoma
