/**
 * @file   decode.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_ARROW_DECODE_H
#define COMMON_ARROW_DECODE_H

#include <string>
#include <unordered_map>

namespace tiledbsoma::common::arrow {
std::unordered_map<std::string, std::string> metadata_string_to_map(const char* metadata_str);
}

#endif