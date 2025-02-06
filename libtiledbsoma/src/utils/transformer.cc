#include "transformer.h"
#include "../soma/soma_coordinates.h"
#include "../geometry/geometry.h"
#include "../geometry/operators/envelope.h"
#include "../geometry/operators/io/write.h"

#include <unordered_set>

namespace tiledbsoma::transformer {
TransformerPipeline::TransformerPipeline(
    std::unique_ptr<ArrowArray> array, std::unique_ptr<ArrowSchema> schema)
    : array(std::move(array))
    , schema(std::move(schema)) {
}

TransformerPipeline::~TransformerPipeline() {
}

ArrowTable TransformerPipeline::asTable() {
    return std::make_pair(std::move(array), std::move(schema));
}

void OutlineTransformer(ArrowArray* array, ArrowSchema* schema, const tiledbsoma::SOMACoordinateSpace& coordinate_space) {
    // Initialize a vector to hold all the Arrow tables containing the
    // transformed geometry data
    std::vector<ArrowTable> tables;
    tables.push_back(ArrowTable(
        std::make_unique<ArrowArray>(ArrowArray{}),
        std::make_unique<ArrowSchema>(ArrowSchema{})));

    NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
        tables.front().first.get(), ArrowType::NANOARROW_TYPE_LARGE_BINARY));
    NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
        tables.front().second.get(), ArrowType::NANOARROW_TYPE_LARGE_BINARY));
    NANOARROW_THROW_NOT_OK(
        ArrowSchemaSetName(tables.front().second.get(), SOMA_GEOMETRY_COLUMN_NAME.c_str()));

    for (size_t i = 0; i < coordinate_space.size(); ++i) {
        const auto axis = coordinate_space.axis(i);

        // Min spatial axis
        tables.push_back(ArrowTable(
            std::make_unique<ArrowArray>(ArrowArray{}),
            std::make_unique<ArrowSchema>(ArrowSchema{})));
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            tables.back().first.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
            tables.back().second.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaSetName(
            tables.back().second.get(),
            (SOMA_GEOMETRY_DIMENSION_PREFIX + axis.name + "__min").c_str()));

        // Max spatial axis
        tables.push_back(ArrowTable(
            std::make_unique<ArrowArray>(), std::make_unique<ArrowSchema>()));
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            tables.back().first.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
            tables.back().second.get(), ArrowType::NANOARROW_TYPE_DOUBLE));
        NANOARROW_THROW_NOT_OK(ArrowSchemaSetName(
            tables.back().second.get(),
            (SOMA_GEOMETRY_DIMENSION_PREFIX + axis.name + "__max").c_str()));
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
        ArrowArrayReserve(tables.front().first.get(), wkb_buffer_size));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayStartAppending(tables.front().first.get()));
    for (size_t i = 1; i < tables.size(); ++i) {
        NANOARROW_THROW_NOT_OK(
            ArrowArrayReserve(tables[i].first.get(), array->length));
        NANOARROW_THROW_NOT_OK(ArrowArrayStartAppending(tables[i].first.get()));
    }

    for (const auto& geometry : geometries) {
        geometry::BinaryBuffer wkb = geometry::to_wkb(geometry);
        geometry::Envelope envelope = geometry::envelope(geometry);

        ArrowBufferView wkb_view;
        wkb_view.data.data = wkb.data();
        wkb_view.size_bytes = static_cast<int64_t>(wkb.size());

        NANOARROW_THROW_NOT_OK(
            ArrowArrayAppendBytes(tables.front().first.get(), wkb_view));

        for (size_t i = 0; i < coordinate_space.size(); ++i) {
            NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(
                tables[2 * i + 1].first.get(), envelope.range.at(i).first));
            NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(
                tables[2 * i + 2].first.get(), envelope.range.at(i).second));
        }
    }

    for (size_t i = 0; i < tables.size(); ++i) {
        ArrowError error;
        NANOARROW_THROW_NOT_OK(
            ArrowArrayFinishBuildingDefault(tables[i].first.get(), &error));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////

    std::unordered_set<std::string> unique_column_names;
    std::unique_ptr<ArrowSchema> arrow_schema = std::make_unique<ArrowSchema>(
        ArrowSchema{});
    std::unique_ptr<ArrowArray> arrow_array = std::make_unique<ArrowArray>(
        ArrowArray{});

    for (int64_t i = 0; i < original_data.second->n_children; ++i) {
        unique_column_names.insert(original_data.second->children[i]->name);
    }
    for (size_t i = 0; i < wkb_data.size(); ++i) {
        unique_column_names.insert(wkb_data[i].second->name);
    }

    NANOARROW_THROW_NOT_OK(ArrowSchemaInitFromType(
        arrow_schema.get(), ArrowType::NANOARROW_TYPE_STRUCT));
    NANOARROW_THROW_NOT_OK(ArrowSchemaAllocateChildren(
        arrow_schema.get(), unique_column_names.size()));
    NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
        arrow_array.get(), ArrowType::NANOARROW_TYPE_STRUCT));
    NANOARROW_THROW_NOT_OK(ArrowArrayAllocateChildren(
        arrow_array.get(), unique_column_names.size()));

    // First add the wkb data columns so that already existing columns in the
    // original data except `soma_geometry` can overwrite the generated columns.

    for (size_t i = 0; i < wkb_data.size(); ++i) {
        ArrowSchemaMove(wkb_data[i].second.get(), arrow_schema->children[i]);
        ArrowArrayMove(wkb_data[i].first.get(), arrow_array->children[i]);
    }

    int64_t index = wkb_data.size();
    for (int64_t i = 0; i < original_data.second->n_children; ++i) {
        if (strcmp(original_data.second->children[i]->name, "soma_geometry") ==
            0) {
            continue;
        }

        bool replaced = false;
        for (size_t j = 0; j < wkb_data.size(); ++j) {
            if (strcmp(
                    arrow_schema->children[j]->name,
                    original_data.second->children[i]->name) == 0) {
                arrow_schema->children[j]->release(arrow_schema->children[j]);
                arrow_array->children[j]->release(arrow_array->children[j]);

                ArrowSchemaMove(
                    original_data.second->children[i],
                    arrow_schema->children[j]);
                ArrowArrayMove(
                    original_data.first->children[i], arrow_array->children[j]);

                replaced = true;
                break;
            }
        }

        if (!replaced) {
            ArrowSchemaMove(
                original_data.second->children[i],
                arrow_schema->children[index]);
            ArrowArrayMove(
                original_data.first->children[i], arrow_array->children[index]);

            ++index;
        }
    }

}

}  // namespace tiledbsoma::transformer