/**
 * @file   utils.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_METADATA_UTILS_H
#define COMMON_METADATA_UTILS_H

#include "../datatype/datatype.h"
#include "metadata.h"

namespace tiledbsoma::common {
/**
 * @brief Decode a TileDB formatted metadata value to a typed value for easier access.
 * 
 * @param type The DataType of the metadata value.
 * @param elements The number of elements of the given type the metadata value consists of. For strings 
 *  this value matches the number of characters excluding the null terminator.
 * @param data A raw pointer the the metadata value data buffer.
 */
MetadataValue decode_metadata(DataType type, uint32_t elements, const void* data);
}  // namespace tiledbsoma::common

#endif