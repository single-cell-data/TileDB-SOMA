#include "soma_geometry_dimension.h"

namespace tiledbsoma {

SOMAGeometryColumn SOMAGeometryColumn::create(
    const Context& ctx, Array& array) {
    std::vector<Dimension> dims;

    const ArraySchema& schema = array.schema();

    for (size_t i = 0; i < schema.domain().ndim(); ++i) {
        if (schema.domain().dimension(i).name().rfind(
                SOMA_GEOMETRY_DIMENSION_PREFIX, 0) == 0) {
            dims.push_back(schema.domain().dimension(i));
        }
    }

    if (dims.size() != 4 && dims.size() != 6) {
        throw TileDBSOMAError(fmt::format(
            "[SOMAGeometryColumn] Spatial dimension invalid count. Expected 4 "
            "or 6, found {}",
            dims.size()));
    }

    for (size_t i = 0; i < schema.attribute_num(); ++i) {
        if (schema.attribute(i).name() == SOMA_GEOMETRY_COLUMN_NAME) {
            return SOMAGeometryColumn(ctx, array, dims, schema.attribute(i));
        }
    }

    throw TileDBSOMAError(fmt::format(
        "[SOMAGeometryColumn] Missing {} attribute",
        SOMA_GEOMETRY_COLUMN_NAME));
}

void SOMAGeometryColumn::_set_dim_ranges(
    ManagedQuery& query, const std::any& ranges) const {
    std::vector<std::pair<double_t, double_t>>
        transformed_ranges = _transform_ranges(
            std::any_cast<std::vector<
                std::pair<std::vector<double_t>, std::vector<double_t>>>>(
                ranges));

    auto domain_limits = _limits();

    // Create a rnage object and reuse if for all dimensions
    std::vector<std::pair<double_t, double_t>> range(1);

    for (size_t i = 0; i < transformed_ranges.size(); ++i) {
        // TODO: Maybe remove?
        impl::type_check<double_t>(dimensions[2 * i].type());
        impl::type_check<double_t>(dimensions[2 * i + 1].type());

        range[0] = std::make_pair(
            domain_limits[i].first,
            std::min(transformed_ranges[i].second, domain_limits[i].second));
        query.select_ranges(dimensions[2 * i].name(), range);

        range[0] = std::make_pair(
            std::max(transformed_ranges[i].first, domain_limits[i].first),
            domain_limits[i].second);
        query.select_ranges(dimensions[2 * i + 1].name(), range);
    }
}

std::vector<std::pair<double_t, double_t>> SOMAGeometryColumn::_limits() const {
    std::vector<std::pair<double_t, double_t>> limits;

    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        if (ArraySchemaExperimental::current_domain(ctx, array.schema())
                .is_empty()) {
            std::pair<double_t, double_t> domain = dimensions[2 * i]
                                                       .domain<double_t>();

            limits.push_back(std::make_pair(domain.first, domain.second));
        } else {
            std::array<double_t, 2>
                domain = ArraySchemaExperimental::current_domain(
                             ctx, array.schema())
                             .ndrectangle()
                             .range<double_t>(dimensions.at(2 * i).name());

            limits.push_back(std::make_pair(domain[0], domain[1]));
        }
    }

    return limits;
}

std::vector<std::pair<double_t, double_t>>
SOMAGeometryColumn::_transform_ranges(
    const std::vector<std::pair<std::vector<double_t>, std::vector<double_t>>>&
        ranges) const {
    if (ranges.size() != 1) {
        throw TileDBSOMAError(
            "Multi ranges are not supported for geometry dimension");
    }

    std::vector<std::pair<double_t, double_t>> transformed_ranges;
    std::vector<double_t> min_ranges = ranges.front().first;
    std::vector<double_t> max_ranges = ranges.front().second;
    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        transformed_ranges.push_back(
            std::make_pair(min_ranges[i], max_ranges[i]));
        transformed_ranges.push_back(
            std::make_pair(min_ranges[i], max_ranges[i]));
    }

    return transformed_ranges;
}

std::any SOMAGeometryColumn::_core_domain_slot() const {
    std::vector<double_t> min, max;
    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        std::pair<double_t, double_t> domain = dimensions[2 * i]
                                                   .domain<double_t>();

        min.push_back(domain.first);
        max.push_back(domain.second);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
};

std::any SOMAGeometryColumn::_non_empty_domain_slot() const {
    std::vector<double_t> min, max;
    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        std::pair<double_t, double_t>
            min_domain = array.non_empty_domain<double_t>(
                dimensions[2 * i].name());
        std::pair<double_t, double_t>
            max_domain = array.non_empty_domain<double_t>(
                dimensions[2 * i + 1].name());

        min.push_back(min_domain.first);
        max.push_back(max_domain.second);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
}

std::any SOMAGeometryColumn::_core_current_domain_slot() const {
    std::vector<double_t> min, max;
    CurrentDomain
        current_domain = tiledb::ArraySchemaExperimental::current_domain(
            ctx, array.schema());
    NDRectangle ndrect = current_domain.ndrectangle();

    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        std::array<double_t, 2> domain = ndrect.range<double_t>(
            dimensions[2 * i].name());

        min.push_back(domain[0]);
        max.push_back(domain[1]);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
}

}  // namespace tiledbsoma