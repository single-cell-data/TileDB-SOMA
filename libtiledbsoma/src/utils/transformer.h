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
 * class provides methods to transform data stored in Arrow Tables. The
 * transformation can be in place to the same Arrow Table or return a new Arrow
 * Table.
 */

#include "arrow_adapter.h"

#include <functional>

#ifndef SOMA_TRANSFORMER_H
#define SOMA_TRANSFORMER_H

namespace tiledbsoma::transformer {
class SOMACoordinateSpace;

class TransformerPipeline {
   public:
    TransformerPipeline(
        std::unique_ptr<ArrowArray> array, std::unique_ptr<ArrowSchema> schema);
    ~TransformerPipeline();

    template <class... Ts>
    TransformerPipeline& transform(
        std::function<void(ArrowArray*, ArrowSchema*, Ts...)> transformer,
        Ts... args) {
        transformer(array.get(), schema.get(), args...);

        return *this;
    }

    ArrowTable asTable();

    std::unique_ptr<ArrowArray> array;
    std::unique_ptr<ArrowSchema> schema;
};

void OutlineTransformer(ArrowArray* array, ArrowSchema* schema, const SOMACoordinateSpace& coordinate_space);

}  // namespace tiledbsoma::transformer

#endif