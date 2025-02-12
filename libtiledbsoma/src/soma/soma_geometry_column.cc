/**
 * @file   soma_geometry_column.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAGeometryColumn class.
 */

#include "soma_geometry_column.h"

namespace tiledbsoma {
std::shared_ptr<SOMAColumn> SOMAGeometryColumn::deserialize(
    const nlohmann::json& soma_schema,
    const Context&,
    const Array& array,
    const std::map<std::string, tiledbsoma::MetadataValue>& metadata) {
    if (!soma_schema.contains(TILEDB_SOMA_SCHEMA_COL_DIM_KEY)) {
        throw TileDBSOMAError(
            "[SOMAGeometryColumn][deserialize] Missing required field "
            "'tiledb_dimensions'");
    }
    if (!soma_schema.contains(TILEDB_SOMA_SCHEMA_COL_ATTR_KEY)) {
        throw TileDBSOMAError(
            "[SOMAGeometryColumn][deserialize] Missing required field "
            "'tiledb_attributes'");
    }

    std::vector<std::string>
        dimension_names = soma_schema[TILEDB_SOMA_SCHEMA_COL_DIM_KEY]
                              .template get<std::vector<std::string>>();

    std::vector<std::string>
        attribute_names = soma_schema[TILEDB_SOMA_SCHEMA_COL_ATTR_KEY]
                              .template get<std::vector<std::string>>();

    if (dimension_names.size() % 2 != 0) {
        throw TileDBSOMAError(std::format(
            "[SOMAGeometryColumn][deserialize] Invalid number of dimensions: "
            "expected number divisible by 2, got {}",
            dimension_names.size()));
    }
    if (attribute_names.size() != 1) {
        throw TileDBSOMAError(std::format(
            "[SOMAGeometryColumn][deserialize] Invalid number of attributes: "
            "expected 1, got {}",
            attribute_names.size()));
    }

    std::vector<Dimension> dimensions;
    std::for_each(
        dimension_names.cbegin(),
        dimension_names.cend(),
        [&array, &dimensions](const std::string& name) {
            dimensions.push_back(array.schema().domain().dimension(name));
        });

    auto attribute = array.schema().attribute(attribute_names[0]);

    if (!metadata.contains(SOMA_COORDINATE_SPACE_KEY)) {
        throw TileDBSOMAError(std::format(
            "[SOMAGeometryColumn][deserialize] Missing required '{}' "
            "metadata key.",
            SOMA_COORDINATE_SPACE_KEY));
    }

    auto coordinate_space = std::apply(
        SOMACoordinateSpace::from_metadata,
        metadata.at(SOMA_COORDINATE_SPACE_KEY));

    return std::make_shared<SOMAGeometryColumn>(
        dimensions, attribute, coordinate_space);
}

std::shared_ptr<SOMAGeometryColumn> SOMAGeometryColumn::create(
    std::shared_ptr<Context> ctx,
    ArrowSchema* schema,
    ArrowSchema* spatial_schema,
    ArrowArray* spatial_array,
    const SOMACoordinateSpace& coordinate_space,
    const std::string& soma_type,
    std::string_view type_metadata,
    PlatformConfig platform_config) {
    std::vector<Dimension> dims;
    if (type_metadata.compare("WKB") != 0) {
        throw TileDBSOMAError(std::format(
            "[SOMAGeometryColumn] "
            "Unkwown type metadata for `{}`: "
            "Expected 'WKB', got '{}'",
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
        SOMAGeometryColumn(dims, attribute.first, coordinate_space));
}

void SOMAGeometryColumn::_set_dim_points(
    ManagedQuery& query, const std::any& points) const {
    std::vector<std::pair<double_t, double_t>>
        transformed_points = _transform_points(
            std::any_cast<std::span<const std::vector<double_t>>>(points));

    // The limits of the current domain if it exists or the core domain
    // otherwise.
    auto limits = _limits(*query.ctx(), *query.schema());

    // Create a range object and reuse if for all dimensions
    std::vector<std::pair<double_t, double_t>> range(1);
    size_t dimensionality = dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS;

    for (size_t i = 0; i < transformed_points.size(); ++i) {
        range[0] = std::make_pair(
            limits[i].first,
            std::min(transformed_points[i].second, limits[i].second));
        query.select_ranges(dimensions[i].name(), range);

        range[0] = std::make_pair(
            std::max(transformed_points[i].first, limits[i].first),
            limits[i].second);
        query.select_ranges(dimensions[i + dimensionality].name(), range);
    }
}

void SOMAGeometryColumn::_set_dim_ranges(
    ManagedQuery& query, const std::any& ranges) const {
    std::vector<std::pair<double_t, double_t>>
        transformed_ranges = _transform_ranges(
            std::any_cast<std::vector<
                std::pair<std::vector<double_t>, std::vector<double_t>>>>(
                ranges));

    // The limits of the current domain if it exists or the core domain
    // otherwise.
    auto limits = _limits(*query.ctx(), *query.schema());

    // Create a range object and reuse if for all dimensions
    std::vector<std::pair<double_t, double_t>> range(1);
    size_t dimensionality = dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS;

    for (size_t i = 0; i < transformed_ranges.size(); ++i) {
        range[0] = std::make_pair(
            limits[i].first,
            std::min(transformed_ranges[i].second, limits[i].second));
        query.select_ranges(dimensions[i].name(), range);

        range[0] = std::make_pair(
            std::max(transformed_ranges[i].first, limits[i].first),
            limits[i].second);
        query.select_ranges(dimensions[i + dimensionality].name(), range);
    }
}

void SOMAGeometryColumn::_set_current_domain_slot(
    NDRectangle& rectangle,
    std::span<const std::any> new_current_domain) const {
    if (TDB_DIM_PER_SPATIAL_AXIS * new_current_domain.size() !=
        dimensions.size()) {
        throw TileDBSOMAError(std::format(
            "[SOMAGeometryColumn] Dimension - Current Domain mismatch. "
            "Expected current domain of size {}, found {}",
            dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS,
            new_current_domain.size()));
    }

    for (size_t i = 0; i < new_current_domain.size(); ++i) {
        auto range = std::any_cast<std::array<double_t, 2>>(
            new_current_domain[i]);
        rectangle.set_range<double_t>(dimensions[i].name(), range[0], range[1]);
    }

    for (size_t i = 0; i < new_current_domain.size(); ++i) {
        auto range = std::any_cast<std::array<double_t, 2>>(
            new_current_domain[i]);
        rectangle.set_range<double_t>(
            dimensions[i + new_current_domain.size()].name(),
            range[0],
            range[1]);
    }
}

std::pair<bool, std::string> SOMAGeometryColumn::_can_set_current_domain_slot(
    std::optional<NDRectangle>& rectangle,
    std::span<const std::any> new_current_domain) const {
    if (new_current_domain.size() !=
        dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS) {
        throw TileDBSOMAError(std::format(
            "[SOMADimension][_can_set_current_domain_slot] Expected current "
            "domain "
            "size is 2, found {}",
            new_current_domain.size()));
    }

    for (size_t i = 0; i < new_current_domain.size(); ++i) {
        auto range = std::any_cast<std::array<double_t, 2>>(
            new_current_domain[i]);

        if (range[0] > range[1]) {
            return std::pair(
                false,
                std::format(
                    "index-column name {}: new lower > new upper",
                    dimensions[i].name()));
        }

        auto dimension_min = dimensions[i];
        auto dimension_max =
            dimensions[i + dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS];

        if (rectangle.has_value()) {
            auto range_min = rectangle.value().range<double_t>(
                dimension_min.name());
            auto range_max = rectangle.value().range<double_t>(
                dimension_max.name());

            if (range[0] > range_min[0]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new lower > old lower (downsize "
                        "is unsupported)",
                        dimension_min.name()));
            }
            if (range[0] > range_max[0]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new lower > old lower (downsize "
                        "is unsupported)",
                        dimension_max.name()));
            }
            if (range[1] < range_min[1]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new upper < old upper (downsize "
                        "is unsupported)",
                        dimension_min.name()));
            }
            if (range[1] < range_max[1]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new upper < old upper (downsize "
                        "is unsupported)",
                        dimension_max.name()));
            }
        } else {
            auto core_domain = std::any_cast<
                std::pair<std::vector<double_t>, std::vector<double_t>>>(
                _core_domain_slot());

            if (range[0] > core_domain.first[i]) {
                return std::pair(
                    false,
                    std::format(
                        "index-column name {}: new lower < limit lower",
                        dimension_min.name()));
            }
            if (range[1] < core_domain.second[i]) {
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
    const Context& ctx, const ArraySchema& schema) const {
    std::vector<std::pair<double_t, double_t>> limits;

    if (ArraySchemaExperimental::current_domain(ctx, schema).is_empty()) {
        for (size_t i = 0; i < dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS;
             ++i) {
            std::pair<double_t, double_t> core_domain = dimensions[i]
                                                            .domain<double_t>();

            limits.push_back(
                std::make_pair(core_domain.first, core_domain.second));
        }
    } else {
        NDRectangle ndrect = ArraySchemaExperimental::current_domain(
                                 ctx, schema)
                                 .ndrectangle();
        for (size_t i = 0; i < dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS;
             ++i) {
            std::array<double_t, 2> range = ndrect.range<double_t>(
                dimensions.at(i).name());

            limits.push_back(std::make_pair(range[0], range[1]));
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
            "Multiranges are not supported for geometry dimension");
    }

    std::vector<std::pair<double_t, double_t>> transformed_ranges;
    std::vector<double_t> min_ranges = ranges.front().first;
    std::vector<double_t> max_ranges = ranges.front().second;
    for (size_t i = 0; i < dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS; ++i) {
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
            "Multipoints are not supported for geometry dimension");
    }

    std::vector<std::pair<double_t, double_t>> transformed_ranges;
    for (size_t i = 0; i < dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS; ++i) {
        transformed_ranges.push_back(
            std::make_pair(points.front()[i], points.front()[i]));
    }

    return transformed_ranges;
}

std::any SOMAGeometryColumn::_core_domain_slot() const {
    std::vector<double_t> min, max;
    for (size_t i = 0; i < dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS; ++i) {
        std::pair<double_t, double_t> core_domain = dimensions[i]
                                                        .domain<double_t>();

        min.push_back(core_domain.first);
        max.push_back(core_domain.second);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
};

std::any SOMAGeometryColumn::_non_empty_domain_slot(Array& array) const {
    std::vector<double_t> min, max;
    size_t dimensionality = dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS;
    for (size_t i = 0; i < dimensionality; ++i) {
        std::pair<double_t, double_t>
            min_non_empty_dom = array.non_empty_domain<double_t>(
                dimensions[i].name());
        std::pair<double_t, double_t>
            max_non_empty_dom = array.non_empty_domain<double_t>(
                dimensions[i + dimensionality].name());

        min.push_back(min_non_empty_dom.first);
        max.push_back(max_non_empty_dom.second);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
}

std::any SOMAGeometryColumn::_non_empty_domain_slot_opt(
    const SOMAContext& ctx, Array& array) const {
    std::vector<double_t> min, max;
    size_t dimensionality = dimensions.size() / 2;
    int32_t is_empty;
    double_t fixed_ned[2];

    for (size_t i = 0; i < dimensionality; ++i) {
        ctx.tiledb_ctx()->handle_error(
            tiledb_array_get_non_empty_domain_from_name(
                ctx.tiledb_ctx()->ptr().get(),
                array.ptr().get(),
                dimensions[i].name().c_str(),  // Min dimension
                fixed_ned,
                &is_empty));

        if (is_empty) {
            return std::make_any<std::optional<
                std::pair<std::vector<double_t>, std::vector<double_t>>>>(
                std::nullopt);
        }

        min.push_back(fixed_ned[0]);

        ctx.tiledb_ctx()->handle_error(
            tiledb_array_get_non_empty_domain_from_name(
                ctx.tiledb_ctx()->ptr().get(),
                array.ptr().get(),
                dimensions[i].name().c_str(),  // Max dimension
                fixed_ned,
                &is_empty));

        if (is_empty) {
            return std::make_any<std::optional<
                std::pair<std::vector<double_t>, std::vector<double_t>>>>(
                std::nullopt);
        }

        min.push_back(fixed_ned[1]);
    }

    return std::make_any<
        std::optional<std::pair<std::vector<double_t>, std::vector<double_t>>>>(
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

    for (size_t i = 0; i < dimensions.size() / TDB_DIM_PER_SPATIAL_AXIS; ++i) {
        std::array<double_t, 2> range = ndrect.range<double_t>(
            dimensions[i].name());

        min.push_back(range[0]);
        max.push_back(range[1]);
    }

    return std::make_any<
        std::pair<std::vector<double_t>, std::vector<double_t>>>(
        std::make_pair(min, max));
}

std::pair<ArrowArray*, ArrowSchema*> SOMAGeometryColumn::arrow_domain_slot(
    const SOMAContext& ctx, Array& array, enum Domainish kind) const {
    switch (domain_type().value()) {
        case TILEDB_FLOAT64: {
            auto parent_schema = ArrowAdapter::make_arrow_schema_parent(
                TDB_DIM_PER_SPATIAL_AXIS, name());
            auto parent_array = ArrowAdapter::make_arrow_array_parent(
                TDB_DIM_PER_SPATIAL_AXIS);

            parent_array->length = 2;
            parent_array->n_buffers = 1;
            parent_array->buffers = (const void**)malloc(sizeof(void*));
            parent_array->buffers[0] = nullptr;

            auto kind_domain = domain_slot<std::vector<double_t>>(
                ctx, array, kind);

            for (size_t i = 0; i < TDB_DIM_PER_SPATIAL_AXIS; ++i) {
                // Generate coordinate space axie schema
                auto child_schema = static_cast<ArrowSchema*>(
                    malloc(sizeof(ArrowSchema)));
                child_schema->format = strdup(
                    ArrowAdapter::to_arrow_format(TILEDB_FLOAT64).data());
                child_schema->name = strdup(
                    coordinate_space.axis(i).name.c_str());
                child_schema->metadata = nullptr;
                child_schema->flags = 0;
                child_schema->n_children = 0;
                child_schema->children = nullptr;
                child_schema->dictionary = nullptr;
                child_schema->release = &ArrowAdapter::release_schema;
                child_schema->private_data = nullptr;

                parent_schema->children[i] = child_schema;

                parent_array->children[i] =
                    ArrowAdapter::make_arrow_array_child<double_t>(std::vector(
                        {kind_domain.first[i], kind_domain.second[i]}));
            }

            return std::make_pair(
                parent_array.release(), parent_schema.release());
        } break;
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
    const SOMAContext& ctx, Array& array) const {
    return ArrowAdapter::arrow_schema_from_tiledb_attribute(
               attribute, *ctx.tiledb_ctx(), array)
        .release();
}

void SOMAGeometryColumn::serialize(nlohmann::json& columns_schema) const {
    nlohmann::json column;

    column[TILEDB_SOMA_SCHEMA_COL_TYPE_KEY] = static_cast<uint32_t>(
        soma_column_datatype_t::SOMA_COLUMN_GEOMETRY);

    column[TILEDB_SOMA_SCHEMA_COL_DIM_KEY] = nlohmann::json::array();
    std::for_each(
        dimensions.cbegin(),
        dimensions.cend(),
        [&column](const Dimension& dim) {
            column[TILEDB_SOMA_SCHEMA_COL_DIM_KEY].push_back(dim.name());
        });
    column[TILEDB_SOMA_SCHEMA_COL_ATTR_KEY] = {attribute.name()};

    columns_schema.push_back(column);
}
}  // namespace tiledbsoma
