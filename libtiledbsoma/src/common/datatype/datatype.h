/**
 * @file   datatype.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_DATATYPE_H
#define COMMON_DATATYPE_H

namespace tiledbsoma::common {
enum class DataType {
    BOOL,
    /** Integral types */
    INT8,
    INT16,
    INT32,
    INT64,
    UINT8,
    UINT16,
    UINT32,
    UINT64,
    /** Floating point types */
    FLOAT32,
    FLOAT64,
    /** Temporal types */
    DATETIME_SEC,
    DATETIME_MS,
    DATETIME_US,
    DATETIME_NS,
    /** String types */
    CHAR,
    STRING_ASCII,
    STRING_UTF8,
    /** Binary types */
    BLOB,
    /** Geometry types */
    GEOM_WKB,
    GEOM_WKT
};
}

#endif