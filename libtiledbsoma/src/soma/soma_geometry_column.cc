#include "soma_geometry_column.h"

namespace tiledbsoma {

std::shared_ptr<SOMAGeometryColumn> SOMAGeometryColumn::create(
    std::shared_ptr<Context> ctx,
    ArrowSchema* schema,
    ArrowSchema* spatial_schema,
    ArrowArray* spatial_array,
    const std::string& soma_type,
    std::string_view type_metadata,
    PlatformConfig platform_config) {
    std::vector<Dimension> dims;
    if (type_metadata.compare("WKB") != 0) {
        throw TileDBSOMAError(std::format(
            "[SOMAGeometryColumn] "
            "Unkwown type metadata for `{}`: "
            "Expected 'WKB', got {}",
            SOMA_GEOMETRY_COLUMN_NAME,
            type_metadata));
    }

    for (int64_t j = 0; j < spatial_schema->n_children; ++j) {
        dims.push_back(ArrowAdapter::tiledb_dimension_from_arrow_schema(
            ctx,
            spatial_schema->children[j],
            spatial_array->children[j],
            soma_type,
            type_metadata,
            SOMA_GEOMETRY_DIMENSION_PREFIX,
            "__min",
            platform_config));
    }

    for (int64_t j = 0; j < spatial_schema->n_children; ++j) {
        dims.push_back(ArrowAdapter::tiledb_dimension_from_arrow_schema(
            ctx,
            spatial_schema->children[j],
            spatial_array->children[j],
            soma_type,
            type_metadata,
            SOMA_GEOMETRY_DIMENSION_PREFIX,
            "__max",
            platform_config));
    }

    auto attribute = ArrowAdapter::tiledb_attribute_from_arrow_schema(
        ctx, schema, type_metadata, platform_config);

    return std::make_shared<SOMAGeometryColumn>(
        SOMAGeometryColumn(dims, attribute.first));
}

void SOMAGeometryColumn::_set_dim_points(
    const std::unique_ptr<ManagedQuery>& query,
    const SOMAContext& ctx,
    const std::any& points) const {
    std::vector<std::pair<double_t, double_t>>
        transformed_points = _transform_points(
            std::any_cast<std::span<const std::vector<double_t>>>(points));

    auto domain_limits = _limits(ctx, *query->schema());

    // Create a range object and reuse if for all dimensions
    std::vector<std::pair<double_t, double_t>> range(1);
    size_t dimensionality = dimensions.size() / 2;

    for (size_t i = 0; i < transformed_points.size(); ++i) {
        range[0] = std::make_pair(
            domain_limits[i].first,
            std::min(transformed_points[i].second, domain_limits[i].second));
        query->select_ranges(dimensions[i].name(), range);

        range[0] = std::make_pair(
            std::max(transformed_points[i].first, domain_limits[i].first),
            domain_limits[i].second);
        query->select_ranges(dimensions[i + dimensionality].name(), range);
    }
}

void SOMAGeometryColumn::_set_dim_ranges(
    const std::unique_ptr<ManagedQuery>& query,
    const SOMAContext& ctx,
    const std::any& ranges) const {
    std::vector<std::pair<double_t, double_t>>
        transformed_ranges = _transform_ranges(
            std::any_cast<std::vector<
                std::pair<std::vector<double_t>, std::vector<double_t>>>>(
                ranges));

    auto domain_limits = _limits(ctx, *query->schema());

    // Create a range object and reuse if for all dimensions
    std::vector<std::pair<double_t, double_t>> range(1);
    size_t dimensionality = dimensions.size() / 2;

    for (size_t i = 0; i < transformed_ranges.size(); ++i) {
        range[0] = std::make_pair(
            domain_limits[i].first,
            std::min(transformed_ranges[i].second, domain_limits[i].second));
        query->select_ranges(dimensions[i].name(), range);

        range[0] = std::make_pair(
            std::max(transformed_ranges[i].first, domain_limits[i].first),
            domain_limits[i].second);
        query->select_ranges(dimensions[i + dimensionality].name(), range);
    }
}

void SOMAGeometryColumn::_set_current_domain_slot(
    NDRectangle& rectangle, std::span<const std::any> domain) const {
    if (2 * domain.size() != dimensions.size()) {
        throw TileDBSOMAError(std::format(
            "[SOMAGeometryColumn] Dimension - Current Domain mismatch. "
            "Expected current domain of size {}, found {}",
            dimensions.size() / 2,
            domain.size()));
    }

    for (size_t i = 0; i < domain.size(); ++i) {
        auto dom = std::any_cast<std::array<double_t, 2>>(domain[i]);
        rectangle.set_range<double_t>(dimensions[i].name(), dom[0], dom[1]);
    }

    for (size_t i = 0; i < domain.size(); ++i) {
        auto dom = std::any_cast<std::array<double_t, 2>>(domain[i]);
        rectangle.set_range<double_t>(
            dimensions[i + domain.size()].name(), dom[0], dom[1]);
    }
}

std::pair<bool, std::string> SOMAGeometryColumn::_can_set_current_domain_slot(
    std::optional<NDRectangle>& rectangle,
    std::span<const std::any> new_domain) const {
    if (new_domain.size() != dimensions.size() / 2) {
        throw TileDBSOMAError(std::format(
            "[SOMADimension][_can_set_current_domain_slot] Expected domain "
            "size is 2, found {}",
            new_domain.size()));
    }

    for (size_t i = 0; i < new_domain.size(); ++i) {
        auto new_dom = std::any_cast<std::array<double_t, 2>>(new_domain[i]);

        if (new_dom[0] > new_dom[1]) {
            return std::pair(
                false,
                std::format(
                    "index-column name {}: new lower > new upper",
                    dimensions[i].name()));
        }

        auto dimension_min = dimensions[i];
        auto dimension_max = dimensions[i + dimensions.size() / 2];

        if (rectangle.has_value()) {
            auto dom_min = rectangle.value().range<double_t>(
                dimension_min.name());
            auto dom_max = rectangle.value().range<double_t>(
                dimension_max.name());

            if (new_dom[0] > dom_min[0]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new lower > old lower (downsize "
                        "is unsupported)",
                        dimension_min.name()));
            }
            if (new_dom[0] > dom_max[0]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new lower > old lower (downsize "
                        "is unsupported)",
                        dimension_max.name()));
            }
            if (new_dom[1] < dom_min[1]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new upper < old upper (downsize "
                        "is unsupported)",
                        dimension_min.name()));
            }
            if (new_dom[1] < dom_max[1]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new upper < old upper (downsize "
                        "is unsupported)",
                        dimension_max.name()));
            }
        } else {
            auto dom = std::any_cast<
                std::pair<std::vector<double_t>, std::vector<double_t>>>(
                _core_domain_slot());

            if (new_dom[0] > dom.first[i]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new lower < limit lower",
                        dimension_min.name()));
            }
            if (new_dom[1] < dom.second[i]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new upper > limit upper",
                        dimension_min.name()));
            }
        }
    }

    return std::pair(true, "");
}

std::vector<std::pair<double_t, double_t>> SOMAGeometryColumn::_limits(
    const SOMAContext& ctx, const ArraySchema& schema) const {
    std::vector<std::pair<double_t, double_t>> limits;

    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        if (ArraySchemaExperimental::current_domain(*ctx.tiledb_ctx(), schema)
                .is_empty()) {
            std::pair<double_t, double_t> domain = dimensions[i]
                                                       .domain<double_t>();

            limits.push_back(std::make_pair(domain.first, domain.second));
        } else {
            std::array<double_t, 2>
                domain = ArraySchemaExperimental::current_domain(
                             *ctx.tiledb_ctx(), schema)
                             .ndrectangle()
                             .range<double_t>(dimensions.at(i).name());

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
    }

    return transformed_ranges;
}

std::vector<std::pair<double_t, double_t>>
SOMAGeometryColumn::_transform_points(
    const std::span<const std::vector<double_t>>& points) const {
    if (points.size() != 1) {
        throw TileDBSOMAError(
            "Multi points are not supported for geometry dimension");
    }

    std::vector<std::pair<double_t, double_t>> transformed_ranges;
    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        transformed_ranges.push_back(
            std::make_pair(points.front()[i], points.front()[i]));
    }

    return transformed_ranges;
}

std::any SOMAGeometryColumn::_core_domain_slot() const {
    std::vector<double_t> min, max;
    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        std::pair<double_t, double_t> domain = dimensions[i].domain<double_t>();

        min.push_back(domain.first);
        max.push_back(domain.second);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
};

std::any SOMAGeometryColumn::_non_empty_domain_slot(Array& array) const {
    std::vector<double_t> min, max;
    size_t dimensionality = dimensions.size() / 2;
    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        std::pair<double_t, double_t>
            min_domain = array.non_empty_domain<double_t>(dimensions[i].name());
        std::pair<double_t, double_t>
            max_domain = array.non_empty_domain<double_t>(
                dimensions[i + dimensionality].name());

        min.push_back(min_domain.first);
        max.push_back(max_domain.second);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
}

std::any SOMAGeometryColumn::_core_current_domain_slot(
    const SOMAContext& ctx, Array& array) const {
    CurrentDomain
        current_domain = tiledb::ArraySchemaExperimental::current_domain(
            *ctx.tiledb_ctx(), array.schema());
    NDRectangle ndrect = current_domain.ndrectangle();

    return _core_current_domain_slot(ndrect);
}

std::any SOMAGeometryColumn::_core_current_domain_slot(
    NDRectangle& ndrect) const {
    std::vector<double_t> min, max;

    for (size_t i = 0; i < dimensions.size() / 2; ++i) {
        std::array<double_t, 2> domain = ndrect.range<double_t>(
            dimensions[i].name());

        min.push_back(domain[0]);
        max.push_back(domain[1]);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
}

ArrowArray* SOMAGeometryColumn::arrow_domain_slot(
    const SOMAContext& ctx, Array& array, enum Domainish kind) const {
    switch (domain_type().value()) {
        case TILEDB_FLOAT64:
            return ArrowAdapter::make_arrow_array_child_var(
                domain_slot<std::vector<double_t>>(ctx, array, kind));
            break;
        default:
            throw TileDBSOMAError(std::format(
                "[SOMAGeometryColumn][arrow_domain_slot] dim {} has unhandled "
                "extended type "
                "{}",
                name(),
                tiledb::impl::type_to_str(domain_type().value())));
    }
}

ArrowSchema* SOMAGeometryColumn::arrow_schema_slot(
    const SOMAContext& ctx, Array& array) {
    return ArrowAdapter::arrow_schema_from_tiledb_attribute(
               attribute, *ctx.tiledb_ctx(), array)
        .release();
}

}  // namespace tiledbsoma