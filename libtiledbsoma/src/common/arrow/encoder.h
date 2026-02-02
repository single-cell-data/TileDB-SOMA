/**
 * @file   encoder.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_ARROW_ENCODER_H
#define COMMON_ARROW_ENCODER_H

#include <string>
#include <unordered_map>

namespace tiledbsoma::common::arrow {
char* metadata_map_to_string(const std::unordered_map<std::string, std::string>& metadata);
}

#endif