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
 * This file defines the Transformer base class. A class extending the
 * Transformer class provides a generic method to transform data stored in Arrow
 * tables. The transformation can be in-place to the same Arrow table, and
 * multiple transformations may be chained with using the TransformerPipeline.
 */

#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include "arrow_adapter.h"

#include <concepts>
#include <functional>
#include <vector>

namespace tiledbsoma {

class Transformer {
   public:
    virtual ~Transformer();

    virtual ArrowTable apply(
        std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>) = 0;
};

class TransformerPipeline {
   public:
    TransformerPipeline(
        std::unique_ptr<ArrowArray> array, std::unique_ptr<ArrowSchema> schema);

    TransformerPipeline(TransformerPipeline&& other);

    virtual ~TransformerPipeline();

    TransformerPipeline& operator=(TransformerPipeline&& other);

    TransformerPipeline& transform(std::shared_ptr<Transformer> transformer) {
        std::tie(array, schema) = transformer->apply(
            std::move(array), std::move(schema));

        return *this;
    }

    template <class T>
        requires std::derived_from<T, Transformer>
    TransformerPipeline& transform(T transformer) {
        std::tie(array, schema) = transformer.apply(
            std::move(array), std::move(schema));

        return *this;
    }

    template <typename T, class... Ts>
        requires std::invocable<
                     T,
                     std::unique_ptr<ArrowArray>,
                     std::unique_ptr<ArrowSchema>,
                     Ts...> &&
                 std::same_as<
                     std::invoke_result_t<
                         T,
                         std::unique_ptr<ArrowArray>,
                         std::unique_ptr<ArrowSchema>,
                         Ts...>,
                     ArrowTable>
    TransformerPipeline& transform(T transformer, Ts... args) {
        std::tie(array, schema) = transformer(
            std::move(array), std::move(schema), args...);

        return *this;
    }

    ArrowTable asTable();

   private:
    std::unique_ptr<ArrowArray> array;
    std::unique_ptr<ArrowSchema> schema;
};
}  // namespace tiledbsoma

#endif