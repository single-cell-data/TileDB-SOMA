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

#include <string_view>

namespace tiledbsoma::common {
using namespace std::string_view_literals;

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
    DATETIME_YEAR,
    DATETIME_MONTH,
    DATETIME_WEEK,
    DATETIME_DAY,
    DATETIME_HR,
    DATETIME_MIN,
    DATETIME_SEC,
    DATETIME_MS,
    DATETIME_US,
    DATETIME_NS,
    DATETIME_PS,
    DATETIME_FS,
    DATETIME_AS,
    TIME_HR,
    TIME_MIN,
    TIME_SEC,
    TIME_MS,
    TIME_US,
    TIME_NS,
    TIME_PS,
    TIME_FS,
    TIME_AS,
    /** String types */
    CHAR,
    STRING_ASCII,
    STRING_UTF8,
    STRING_UTF16,
    STRING_UTF32,
    STRING_UCS2,
    STRING_UCS4,
    /** Binary types */
    BLOB,
    /** Geometry types */
    GEOM_WKB,
    GEOM_WKT
};

constexpr std::string_view getName(DataType type) {
    switch (type) {
        case DataType::BOOL:
            return "BOOL"sv;
        case DataType::INT8:
            return "INT8"sv;
        case DataType::INT16:
            return "INT16"sv;
        case DataType::INT32:
            return "INT32"sv;
        case DataType::INT64:
            return "INT64"sv;
        case DataType::UINT8:
            return "UINT8"sv;
        case DataType::UINT16:
            return "UINT16"sv;
        case DataType::UINT32:
            return "UINT32"sv;
        case DataType::UINT64:
            return "UINT64"sv;
        case DataType::FLOAT32:
            return "FLOAT32"sv;
        case DataType::FLOAT64:
            return "FLOAT64"sv;
        case DataType::DATETIME_YEAR:
            return "DATETIME_YEAR"sv;
        case DataType::DATETIME_MONTH:
            return "DATETIME_MONTH"sv;
        case DataType::DATETIME_WEEK:
            return "DATETIME_WEEK"sv;
        case DataType::DATETIME_DAY:
            return "DATETIME_DAY"sv;
        case DataType::DATETIME_HR:
            return "DATETIME_HR"sv;
        case DataType::DATETIME_MIN:
            return "DATETIME_MIN"sv;
        case DataType::DATETIME_SEC:
            return "DATETIME_SEC"sv;
        case DataType::DATETIME_MS:
            return "DATETIME_MS"sv;
        case DataType::DATETIME_US:
            return "DATETIME_US"sv;
        case DataType::DATETIME_NS:
            return "DATETIME_NS"sv;
        case DataType::DATETIME_PS:
            return "DATETIME_PS"sv;
        case DataType::DATETIME_FS:
            return "DATETIME_FS"sv;
        case DataType::DATETIME_AS:
            return "DATETIME_AS"sv;
        case DataType::TIME_HR:
            return "TIME_HR"sv;
        case DataType::TIME_MIN:
            return "TIME_MIN"sv;
        case DataType::TIME_SEC:
            return "TIME_SEC"sv;
        case DataType::TIME_MS:
            return "TIME_MS"sv;
        case DataType::TIME_US:
            return "TIME_US"sv;
        case DataType::TIME_NS:
            return "TIME_NS"sv;
        case DataType::TIME_PS:
            return "TIME_PS"sv;
        case DataType::TIME_FS:
            return "TIME_FS"sv;
        case DataType::TIME_AS:
            return "TIME_AS"sv;
        case DataType::CHAR:
            return "CHAR"sv;
        case DataType::STRING_ASCII:
            return "STRING_ASCII"sv;
        case DataType::STRING_UTF8:
            return "STRING_UTF8"sv;
        case DataType::STRING_UTF16:
            return "STRING_UTF16"sv;
        case DataType::STRING_UTF32:
            return "STRING_UTF32"sv;
        case DataType::STRING_UCS2:
            return "STRING_UCS2"sv;
        case DataType::STRING_UCS4:
            return "STRING_UCS4"sv;
        case DataType::BLOB:
            return "BLOB"sv;
        case DataType::GEOM_WKB:
            return "GEOM_WKB"sv;
        case DataType::GEOM_WKT:
            return "GEOM_WKT"sv;
    }
}
}  // namespace tiledbsoma::common

#endif