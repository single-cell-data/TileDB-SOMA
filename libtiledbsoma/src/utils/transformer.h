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
#include "common/arrow/utils.h"

#include <concepts>
#include <functional>
#include <vector>

namespace tiledbsoma {

class Transformer {
   public:
    virtual ~Transformer();

    virtual common::arrow::ArrowTable apply(
        common::arrow::managed_unique_ptr<ArrowArray>, common::arrow::managed_unique_ptr<ArrowSchema>) = 0;
};

class TransformerPipeline {
   public:
    TransformerPipeline(
        common::arrow::managed_unique_ptr<ArrowArray> array, common::arrow::managed_unique_ptr<ArrowSchema> schema);

    TransformerPipeline(TransformerPipeline&& other);

    virtual ~TransformerPipeline();

    TransformerPipeline& operator=(TransformerPipeline&& other);

    TransformerPipeline& transform(std::shared_ptr<Transformer> transformer) {
        std::tie(array, schema) = transformer->apply(std::move(array), std::move(schema));

        return *this;
    }

    template <class T>
        requires std::derived_from<T, Transformer>
    TransformerPipeline& transform(T transformer) {
        std::tie(array, schema) = transformer.apply(std::move(array), std::move(schema));

        return *this;
    }

    template <typename T, class... Ts>
        requires std::invocable<
                     T,
                     common::arrow::managed_unique_ptr<ArrowArray>,
                     common::arrow::managed_unique_ptr<ArrowSchema>,
                     Ts...> &&
                 std::same_as<
                     std::invoke_result_t<
                         T,
                         common::arrow::managed_unique_ptr<ArrowArray>,
                         common::arrow::managed_unique_ptr<ArrowSchema>,
                         Ts...>,
                     common::arrow::ArrowTable>
    TransformerPipeline& transform(T transformer, Ts... args) {
        std::tie(array, schema) = transformer(std::move(array), std::move(schema), args...);

        return *this;
    }

    common::arrow::ArrowTable asTable();

   private:
    common::arrow::managed_unique_ptr<ArrowArray> array;
    common::arrow::managed_unique_ptr<ArrowSchema> schema;
};
}  // namespace tiledbsoma

#endif