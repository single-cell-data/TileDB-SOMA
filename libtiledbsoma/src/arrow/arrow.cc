#include "arrow.h"

#include <sparrow/sparrow.hpp>

namespace tiledbsoma::arrow {

std::shared_ptr<const Array> Array::create(
    ArrowSchema* schema, ArrowArray* array) {
    return make_shared_array(schema, array);
}

Array::Array(ArrowSchema* schema, ArrowArray* array)
    : array_(std::make_unique<sparrow::array>(array, schema)) {
}

std::string_view Array::name() const noexcept {
    return array_->name().value_or("");
}

tiledb_datatype_t Array::tdb_type() const noexcept {
    switch (array_->data_type()) {
        case sparrow::data_type::BOOL:
            return TILEDB_BOOL;
        case sparrow::data_type::INT8:
            return TILEDB_INT8;
        case sparrow::data_type::INT16:
            return TILEDB_INT16;
        case sparrow::data_type::INT32:
            return TILEDB_INT32;
        case sparrow::data_type::INT64:
            return TILEDB_INT64;
        case sparrow::data_type::UINT8:
            return TILEDB_UINT8;
        case sparrow::data_type::UINT16:
            return TILEDB_UINT16;
        case sparrow::data_type::UINT32:
            return TILEDB_UINT32;
        case sparrow::data_type::UINT64:
            return TILEDB_UINT64;
        case sparrow::data_type::FLOAT:
            return TILEDB_FLOAT32;
        case sparrow::data_type::DOUBLE:
            return TILEDB_FLOAT64;
        case sparrow::data_type::STRING:
        case sparrow::data_type::LARGE_STRING:
            return TILEDB_STRING_ASCII;
        case sparrow::data_type::BINARY:
        case sparrow::data_type::LARGE_BINARY:
            return TILEDB_BLOB;
        default:
            return TILEDB_ANY;
    }
}

std::pair<ArrowArray, ArrowSchema> Array::export_to_c() const noexcept {
    return sparrow::extract_arrow_structures(std::move(*array_));
}

void Array::export_to_c(ArrowArray* array, ArrowSchema* schema) const noexcept {
    auto [array_temp, schema_temp] = export_to_c();

    *array = array_temp;
    *schema = schema_temp;
}
}  // namespace tiledbsoma::arrow