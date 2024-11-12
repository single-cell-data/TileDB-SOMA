#include "soma_attribute.h"

namespace tiledbsoma {
std::shared_ptr<SOMAAttribute> SOMAAttribute::create(
    std::shared_ptr<Context> ctx,
    ArrowSchema* schema,
    std::string_view type_metadata,
    PlatformConfig platform_config) {
    auto attribute = ArrowAdapter::tiledb_attribute_from_arrow_schema(
        ctx, schema, type_metadata, platform_config);

    auto result = std::make_shared<SOMAAttribute>(
        SOMAAttribute(*ctx, attribute.first));

    result->enumeration = attribute.second;
    return result;
}

void SOMAAttribute::_set_dim_ranges(
    const std::unique_ptr<ManagedQuery>& query, const std::any& ranges) const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}

void SOMAAttribute::_set_current_domain_slot(
    NDRectangle& rectangle, const std::vector<const void*>& domain) const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}

std::any SOMAAttribute::_core_domain_slot() const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}

std::any SOMAAttribute::_non_empty_domain_slot(Array& array) const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}

std::any SOMAAttribute::_core_current_domain_slot(Array& array) const {
    throw TileDBSOMAError(fmt::format(
        "[SOMAAttribute] Column with name {} is not an index column", name()));
}
}  // namespace tiledbsoma
