/**
 * @file   common.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#include "common.h"

#include <cxxabi.h>

namespace tiledbsoma::common {
#pragma region Status implementation

Status::Status()
    : code_(StatusCode::OK) {
}

Status::Status(StatusCode code, std::string_view origin, std::string_view message)
    : code_(code)
    , origin_(origin)
    , message_(message) {
}

size_t get_max_capacity(tiledb_datatype_t index_type) {
    switch (index_type) {
        case TILEDB_INT8:
            return std::numeric_limits<int8_t>::max();
        case TILEDB_UINT8:
            return std::numeric_limits<uint8_t>::max();
        case TILEDB_INT16:
            return std::numeric_limits<int16_t>::max();
        case TILEDB_UINT16:
            return std::numeric_limits<uint16_t>::max();
        case TILEDB_INT32:
            return std::numeric_limits<int32_t>::max();
        case TILEDB_UINT32:
            return std::numeric_limits<uint32_t>::max();
        case TILEDB_INT64:
            return std::numeric_limits<int64_t>::max();
        case TILEDB_UINT64:
            return std::numeric_limits<uint64_t>::max();
        default:
            throw std::runtime_error(
                "[get_max_capacity] Saw invalid enumeration index type when trying to extend enumeration");
    }
}

std::string demangle_name(std::string_view name) {
    int status = 0;
    std::unique_ptr<char, void (*)(void*)> demangled_name{
        abi::__cxa_demangle(name.data(), nullptr, nullptr, &status), std::free};

    if (status != 0) {
        return name.data();
    }

    return demangled_name.get();
}

#pragma endregion
}  // namespace tiledbsoma::common