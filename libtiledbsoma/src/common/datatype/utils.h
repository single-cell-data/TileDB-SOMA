/**
 * @file   utils.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_DATATYPE_UTILS_H
#define COMMON_DATATYPE_UTILS_H

#include "datatype.h"

#include <tiledb/tiledb>

namespace tiledbsoma::common {
DataType tiledb_to_datatype(tiledb_datatype_t type);
}

#endif