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
    boolean,
    /** Integral types */
    int8,
    int16,
    int32,
    int64,
    uint8,
    uint16,
    uint32,
    uint64,
    /** Floating point types */
    float32,
    float64,
    /** Temporal types */
    datetime_year,
    datetime_month,
    datetime_week,
    datetime_day,
    datetime_hr,
    datetime_min,
    datetime_sec,
    datetime_ms,
    datetime_us,
    datetime_ns,
    datetime_ps,
    datetime_fs,
    datetime_as,
    time_hr,
    time_min,
    time_sec,
    time_ms,
    time_us,
    time_ns,
    time_ps,
    time_fs,
    time_as,
    /** String types */
    character,
    string_ascii,
    string_utf8,
    string_utf16,
    string_utf32,
    string_ucs2,
    string_ucs4,
    /** Binary types */
    blob,
    /** Geometry types */
    geom_wkb,
    geom_wkt
};

constexpr std::string_view getName(DataType type) {
    switch (type) {
        case DataType::boolean:
            return "BOOL"sv;
        case DataType::int8:
            return "INT8"sv;
        case DataType::int16:
            return "INT16"sv;
        case DataType::int32:
            return "INT32"sv;
        case DataType::int64:
            return "INT64"sv;
        case DataType::uint8:
            return "UINT8"sv;
        case DataType::uint16:
            return "UINT16"sv;
        case DataType::uint32:
            return "UINT32"sv;
        case DataType::uint64:
            return "UINT64"sv;
        case DataType::float32:
            return "FLOAT32"sv;
        case DataType::float64:
            return "FLOAT64"sv;
        case DataType::datetime_year:
            return "DATETIME_YEAR"sv;
        case DataType::datetime_month:
            return "DATETIME_MONTH"sv;
        case DataType::datetime_week:
            return "DATETIME_WEEK"sv;
        case DataType::datetime_day:
            return "DATETIME_DAY"sv;
        case DataType::datetime_hr:
            return "DATETIME_HR"sv;
        case DataType::datetime_min:
            return "DATETIME_MIN"sv;
        case DataType::datetime_sec:
            return "DATETIME_SEC"sv;
        case DataType::datetime_ms:
            return "DATETIME_MS"sv;
        case DataType::datetime_us:
            return "DATETIME_US"sv;
        case DataType::datetime_ns:
            return "DATETIME_NS"sv;
        case DataType::datetime_ps:
            return "DATETIME_PS"sv;
        case DataType::datetime_fs:
            return "DATETIME_FS"sv;
        case DataType::datetime_as:
            return "DATETIME_AS"sv;
        case DataType::time_hr:
            return "TIME_HR"sv;
        case DataType::time_min:
            return "TIME_MIN"sv;
        case DataType::time_sec:
            return "TIME_SEC"sv;
        case DataType::time_ms:
            return "TIME_MS"sv;
        case DataType::time_us:
            return "TIME_US"sv;
        case DataType::time_ns:
            return "TIME_NS"sv;
        case DataType::time_ps:
            return "TIME_PS"sv;
        case DataType::time_fs:
            return "TIME_FS"sv;
        case DataType::time_as:
            return "TIME_AS"sv;
        case DataType::character:
            return "CHAR"sv;
        case DataType::string_ascii:
            return "STRING_ASCII"sv;
        case DataType::string_utf8:
            return "STRING_UTF8"sv;
        case DataType::string_utf16:
            return "STRING_UTF16"sv;
        case DataType::string_utf32:
            return "STRING_UTF32"sv;
        case DataType::string_ucs2:
            return "STRING_UCS2"sv;
        case DataType::string_ucs4:
            return "STRING_UCS4"sv;
        case DataType::blob:
            return "BLOB"sv;
        case DataType::geom_wkb:
            return "GEOM_WKB"sv;
        case DataType::geom_wkt:
            return "GEOM_WKT"sv;
    }
}
}  // namespace tiledbsoma::common

#endif