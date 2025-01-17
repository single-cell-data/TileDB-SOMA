/**
 * @file   util.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines utilities.
 */

#include "utils/util.h"
#include <cstring>
#include "logger.h"

namespace tiledbsoma::util {

template <typename T>
VarlenBufferPair to_varlen_buffers(std::vector<T> data, bool arrow) {
    size_t nbytes = 0;
    for (auto& elem : data) {
        nbytes += elem.size();
    }

    std::string result;
    std::vector<uint64_t> offsets(data.size() + 1);
    size_t offset = 0;
    size_t idx = 0;

    for (auto& elem : data) {
        result += elem;
        offsets[idx++] = offset;
        offset += elem.size();
    }
    offsets[idx] = offset;

    // Remove extra arrow offset when creating buffers for TileDB write
    if (!arrow) {
        offsets.pop_back();
    }

    return {result, offsets};
}

template VarlenBufferPair to_varlen_buffers(
    std::vector<std::string>, bool arrow);

bool is_tiledb_uri(std::string_view uri) {
    return uri.find("tiledb://") == 0;
}

std::string rstrip_uri(std::string_view uri) {
    return std::regex_replace(std::string(uri), std::regex("/+$"), "");
}

std::vector<uint8_t> cast_bit_to_uint8(ArrowSchema* schema, ArrowArray* array) {
    if (strcmp(schema->format, "b") != 0) {
        throw TileDBSOMAError(std::format(
            "_cast_bit_to_uint8 expected column format to be 'b' but saw {}",
            schema->format));
    }

    uint8_t* data;
    if (array->n_buffers == 3) {
        data = (uint8_t*)array->buffers[2];
    } else {
        data = (uint8_t*)array->buffers[1];
    }

    std::vector<uint8_t> casted(array->length);
    ArrowBitsUnpackInt8(
        data,
        array->offset,
        array->length,
        reinterpret_cast<int8_t*>(casted.data()));
    return casted;
}

std::shared_ptr<SOMAColumn> find_column_by_name(
    std::span<const std::shared_ptr<SOMAColumn>> columns,
    std::string_view name) {
    auto column_it = std::find_if(
        columns.begin(), columns.end(), [&](auto col) {
            return col->name() == name;
        });

    if (column_it == columns.end()) {
        throw TileDBSOMAError(std::format(
            "[ArrowAdapter][tiledb_schema_from_arrow_schema] Index column "
            "'{}' missing",
            name));
    }

    return *column_it;
}

};  // namespace tiledbsoma::util
