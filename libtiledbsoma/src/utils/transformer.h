/**
 * @file   transformer.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the Transformer class. A class implementing the Transformer
 * class provides a generic method to transform data stored in Arrow Tables. The
 * transformation can be in place to the same Arrow Table and multiple
 * transformations can be chained with using the TransformerPipeline.
 */

#include "arrow_adapter.h"

#include <concepts>
#include <functional>
#include <vector>

#ifndef SOMA_TRANSFORMER_H
#define SOMA_TRANSFORMER_H

namespace tiledbsoma::transformer {
class tiledbsoma::SOMACoordinateSpace;

template <class... Ts>
class Transformer {
   public:
    virtual void apply(ArrowArray*, ArrowSchema*, Ts...) = 0;
};

class TransformerPipeline {
   public:
    TransformerPipeline(
        std::unique_ptr<ArrowArray> array, std::unique_ptr<ArrowSchema> schema);
    ~TransformerPipeline();

    template <class T, class... Ts>
        requires std::derived_from<T, Transformer<Ts...>>
    TransformerPipeline& transform(T transformer, Ts... args) {
        transformer.apply(array.get(), schema.get(), args...);

        return *this;
    }

    ArrowTable asTable();

   private:
    std::unique_ptr<ArrowArray> array;
    std::unique_ptr<ArrowSchema> schema;
};

class OutlineTransformer : public Transformer<tiledbsoma::SOMACoordinateSpace> {
   public:
    void apply(
        ArrowArray*, ArrowSchema*, tiledbsoma::SOMACoordinateSpace) override;

   private:
    /**
     * @brief Cast an array containing the outer rings of polygons to an Arrow
     * array holding the WKB encoded polygons and generate the additional index
     * column arrays based on the spatial axes.
     */
    std::vector<ArrowTable> _cast_polygon_vertex_list_to_wkb(
        ArrowArray* array,
        const tiledbsoma::SOMACoordinateSpace& coordinate_space);
};

}  // namespace tiledbsoma::transformer

#endif