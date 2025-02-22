#include "soma_transformers.h"
#include "../geometry/geometry.h"
#include "../geometry/operators/envelope.h"
#include "../geometry/operators/io/write.h"

namespace tiledbsoma {
OutlineTransformer::OutlineTransformer(SOMACoordinateSpace coordinate_space)
    : coordinate_space(coordinate_space) {
}

OutlineTransformer::~OutlineTransformer() {
}

ArrowTable OutlineTransformer::apply(
    std::unique_ptr<ArrowArray> array, std::unique_ptr<ArrowSchema> schema) {
    std::vector<std::unique_ptr<ArrowArray>> generated_arrays;
    std::vector<std::unique_ptr<ArrowSchema>> generated_schemas;

    for (int64_t i = 0; i < schema->n_children; ++i) {
        /**
         * If `soma_geometry` conforms to specific formats, automatically
         * convert to WKB and create additional index columns for spatial axes.
         *
         * If the `soma_geometry` array is a WKB binary, users are expected to
         * provide the additional index columns for spatial axes.
         */

        if (strcmp(schema->children[i]->name, "soma_geometry") == 0 &&
            strcmp(schema->children[i]->format, "+l") == 0) {
            std::tie(generated_arrays, generated_schemas) =
                _cast_polygon_vertex_list_to_wkb(
                    array->children[i], coordinate_space);

            break;
        }
    }

    int64_t soma_gometry_index = -1;
    for (int64_t i = 0; i < schema->n_children; ++i) {
        if (strcmp(
                schema->children[i]->name, SOMA_GEOMETRY_COLUMN_NAME.c_str()) ==
            0) {
            soma_gometry_index = i;
            break;
        }
    }

    if (soma_gometry_index == -1) {
        throw std::runtime_error(std::format(
            "[OutlineTransformer][apply] Missing schema child with name {}",
            SOMA_GEOMETRY_COLUMN_NAME));
    }

    array = ArrowAdapter::arrow_array_remove_at_index(
        std::move(array), soma_gometry_index);
    schema = ArrowAdapter::arrow_schema_remove_at_index(
        std::move(schema), soma_gometry_index);

    array = ArrowAdapter::arrow_array_insert_at_index(
        std::move(array), std::move(generated_arrays), soma_gometry_index);

    schema = ArrowAdapter::arrow_schema_insert_at_index(
        std::move(schema), std::move(generated_schemas), soma_gometry_index);

    return std::make_pair(std::move(array), std::move(schema));
}

std::pair<
    std::vector<std::unique_ptr<ArrowArray>>,
    std::vector<std::unique_ptr<ArrowSchema>>>
OutlineTransformer::_cast_polygon_vertex_list_to_wkb(
    ArrowArray* array,
    const tiledbsoma::SOMACoordinateSpace& coordinate_space) {
    // Initialize a vector to hold all the Arrow tables containing the
    // transformed geometry data
    std::vector<std::unique_ptr<ArrowArray>> arrays;
    std::vector<std::unique_ptr<ArrowSchema>> schemas;

    arrays.push_back(std::make_unique<ArrowArray>(ArrowArray{}));
    schemas.push_back(std::make_unique<ArrowSchema>(ArrowSchema{}));

    NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
        arrays.front().get(), ArrowType::NANOARROW_TYPE_LARGE_BINARY));
    NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
        schemas.front().get(), ArrowType::NANOARROW_TYPE_LARGE_BINARY));
    schemas.front()->name = strdup("soma_geometry");

    for (size_t i = 0; i < coordinate_space.size(); ++i) {
        const auto axis = coordinate_space.axis(i);

        // Min spatial axis
        arrays.push_back(std::move(std::make_unique<ArrowArray>(ArrowArray{})));
        schemas.push_back(
            std::move(std::make_unique<ArrowSchema>(ArrowSchema{})));
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            arrays.back().get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
            schemas.back().get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        schemas.back()->name = strdup(
            (SOMA_GEOMETRY_DIMENSION_PREFIX + axis.name + "__min").c_str());

        // Max spatial axis
        arrays.push_back(std::make_unique<ArrowArray>(ArrowArray{}));
        schemas.push_back(std::make_unique<ArrowSchema>(ArrowSchema{}));
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            arrays.back().get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
            schemas.back().get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        schemas.back()->name = strdup(
            (SOMA_GEOMETRY_DIMENSION_PREFIX + axis.name + "__max").c_str());
    }

    // Large list of doubles
    const uint32_t* offset = static_cast<const uint32_t*>(array->buffers[1]);
    const double_t* data = static_cast<const double_t*>(
        array->children[0]->buffers[1]);

    size_t wkb_buffer_size = 0;
    std::vector<geometry::GenericGeometry> geometries;

    for (int64_t index = 0; index < array->length; ++index) {
        int64_t stop_index = index < array->length - 1 ?
                                 offset[index + 1] :
                                 array->children[0]->length;

        std::vector<geometry::BasePoint> ring;
        for (int64_t j = offset[index]; j < stop_index; j += 2) {
            ring.push_back(geometry::BasePoint(data[j], data[j + 1]));
        }

        geometries.push_back(
            geometry::GenericGeometry(geometry::Polygon(std::move(ring))));
        wkb_buffer_size += wkb_size(geometries.back());
    }

    NANOARROW_THROW_NOT_OK(
        ArrowArrayReserve(arrays.front().get(), wkb_buffer_size));
    NANOARROW_THROW_NOT_OK(ArrowArrayStartAppending(arrays.front().get()));
    for (size_t i = 1; i < arrays.size(); ++i) {
        NANOARROW_THROW_NOT_OK(
            ArrowArrayReserve(arrays[i].get(), array->length));
        NANOARROW_THROW_NOT_OK(ArrowArrayStartAppending(arrays[i].get()));
    }

    for (const auto& geometry : geometries) {
        geometry::BinaryBuffer wkb = geometry::to_wkb(geometry);
        geometry::Envelope envelope = geometry::envelope(geometry);

        ArrowBufferView wkb_view;
        wkb_view.data.data = wkb.data();
        wkb_view.size_bytes = static_cast<int64_t>(wkb.size());

        NANOARROW_THROW_NOT_OK(
            ArrowArrayAppendBytes(arrays.front().get(), wkb_view));

        for (size_t i = 0; i < coordinate_space.size(); ++i) {
            NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(
                arrays[2 * i + 1].get(), envelope.range.at(i).first));
            NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(
                arrays[2 * i + 2].get(), envelope.range.at(i).second));
        }
    }

    for (size_t i = 0; i < arrays.size(); ++i) {
        ArrowError error;
        NANOARROW_THROW_NOT_OK(
            ArrowArrayFinishBuildingDefault(arrays[i].get(), &error));
    }

    return std::make_pair(std::move(arrays), std::move(schemas));
}
}  // namespace tiledbsoma