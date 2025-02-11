#include "transformer.h"

namespace tiledbsoma {
TransformerPipeline::TransformerPipeline(
    std::unique_ptr<ArrowArray> array, std::unique_ptr<ArrowSchema> schema)
    : array(std::move(array))
    , schema(std::move(schema)) {
}

TransformerPipeline::TransformerPipeline(TransformerPipeline&& other)
    : array(std::move(other.array))
    , schema(std::move(other.schema)) {
}

TransformerPipeline::~TransformerPipeline() {
}

TransformerPipeline& TransformerPipeline::operator=(
    TransformerPipeline&& other) {
    if (this != &other) {
        this->array = std::move(other.array);
        this->schema = std::move(other.schema);
    }

    return *this;
}

ArrowTable TransformerPipeline::asTable() {
    return std::make_pair(std::move(array), std::move(schema));
}

Transformer::~Transformer() {
}

}  // namespace tiledbsoma