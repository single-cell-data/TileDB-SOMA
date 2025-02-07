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
    void apply(
        ArrowArray*, ArrowSchema*, tiledbsoma::SOMACoordinateSpace) override;
};

}  // namespace tiledbsoma::transformer

#endif