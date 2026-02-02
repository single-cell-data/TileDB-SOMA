/**
 * @file   encoder.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#include "encoder.h"

#include <cstring>
#include <sstream>

namespace tiledbsoma::common::arrow {
char* metadata_map_to_string(const std::unordered_map<std::string, std::string>& metadata) {
    if (metadata.empty()) {
        return nullptr;
    }

    std::stringstream buffer(std::ios::out | std::ios::binary);

    int32_t metadata_count = metadata.size();
    buffer.write(reinterpret_cast<char*>(std::addressof(metadata_count)), sizeof(int32_t));

    for (const auto& [key, value] : metadata) {
        int32_t key_length = key.size();
        int32_t value_length = value.size();

        buffer.write(reinterpret_cast<char*>(std::addressof(key_length)), sizeof(int32_t));
        buffer.write(key.data(), key_length);
        buffer.write(reinterpret_cast<char*>(std::addressof(value_length)), sizeof(int32_t));
        buffer.write(value.data(), value_length);
    }

    buffer.flush();
    char* result = new char[buffer.view().size()];
    std::memcpy(result, buffer.view().data(), buffer.view().size());

    return result;
}
}  // namespace tiledbsoma::common::arrow