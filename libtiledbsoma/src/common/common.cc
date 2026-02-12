/**
 * @file   common.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#include "common.h"

#if defined(__GNUC__) && !defined(__clang__)
#include <cxxabi.h>
#endif

namespace tiledbsoma::common {

size_t enumeration_value_count(const tiledb::Context& ctx, const tiledb::Enumeration& enumeration) {
    const void* data_buffer;
    uint64_t data_size;
    ctx.handle_error(tiledb_enumeration_get_data(ctx.ptr().get(), enumeration.ptr().get(), &data_buffer, &data_size));

    if (enumeration.cell_val_num() == TILEDB_VAR_NUM) {
        const void* offsets_buffer;
        uint64_t offsets_size;
        ctx.handle_error(
            tiledb_enumeration_get_offsets(ctx.ptr().get(), enumeration.ptr().get(), &offsets_buffer, &offsets_size));

        return offsets_size / sizeof(uint64_t);
    } else {
        return data_size / (tiledb::impl::type_size(enumeration.type()) * enumeration.cell_val_num());
    }
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

#if defined(__GNUC__) && !defined(__clang__)
std::string demangle_name(std::string_view name) {
    int status = 0;
    std::unique_ptr<char, void (*)(void*)> demangled_name{
        abi::__cxa_demangle(name.data(), nullptr, nullptr, &status), std::free};

    if (status != 0) {
        return name.data();
    }

    return demangled_name.get();
}
#else
std::string demangle_name(std::string_view name) {
    return std::string(name);
}
#endif

#pragma endregion
}  // namespace tiledbsoma::common