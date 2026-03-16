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
MetadataValue decode_metadata(DataType type, uint32_t elements, const void* data);
}

#endif