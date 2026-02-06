#include "transformer.h"

namespace tiledbsoma {
TransformerPipeline::TransformerPipeline(
    common::arrow::managed_unique_ptr<ArrowArray> array, common::arrow::managed_unique_ptr<ArrowSchema> schema)
    : array(std::move(array))
    , schema(std::move(schema)) {
}

TransformerPipeline::TransformerPipeline(TransformerPipeline&& other)
    : array(std::move(other.array))
    , schema(std::move(other.schema)) {
}

TransformerPipeline::~TransformerPipeline() {
}

TransformerPipeline& TransformerPipeline::operator=(TransformerPipeline&& other) {
    if (this != &other) {
        this->array = std::move(other.array);
        this->schema = std::move(other.schema);
    }

    return *this;
}

common::arrow::ArrowTable TransformerPipeline::asTable() {
    return std::make_pair(std::move(array), std::move(schema));
}

Transformer::~Transformer() {
}

}  // namespace tiledbsoma