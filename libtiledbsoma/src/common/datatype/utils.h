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

#include <optional>
#include <string_view>
#include <tuple>

namespace tiledbsoma::common::type {
namespace impl {
using namespace std::string_view_literals;
constexpr DataType tiledb_to_datatype(tiledb_datatype_t type) {
    switch (type) {
        case TILEDB_BOOL:
            return DataType::boolean;
        case TILEDB_INT8:
            return DataType::int8;
        case TILEDB_INT16:
            return DataType::int16;
        case TILEDB_INT32:
            return DataType::int32;
        case TILEDB_INT64:
            return DataType::int64;
        case TILEDB_UINT8:
            return DataType::uint8;
        case TILEDB_UINT16:
            return DataType::uint16;
        case TILEDB_UINT32:
            return DataType::uint32;
        case TILEDB_UINT64:
            return DataType::uint64;
        case TILEDB_FLOAT32:
            return DataType::float32;
        case TILEDB_FLOAT64:
            return DataType::float64;
        case TILEDB_DATETIME_YEAR:
            return DataType::datetime_year;
        case TILEDB_DATETIME_MONTH:
            return DataType::datetime_month;
        case TILEDB_DATETIME_WEEK:
            return DataType::datetime_week;
        case TILEDB_DATETIME_DAY:
            return DataType::datetime_day;
        case TILEDB_DATETIME_HR:
            return DataType::datetime_hr;
        case TILEDB_DATETIME_MIN:
            return DataType::datetime_min;
        case TILEDB_DATETIME_SEC:
            return DataType::datetime_sec;
        case TILEDB_DATETIME_MS:
            return DataType::datetime_ms;
        case TILEDB_DATETIME_US:
            return DataType::datetime_us;
        case TILEDB_DATETIME_NS:
            return DataType::datetime_ns;
        case TILEDB_DATETIME_PS:
            return DataType::datetime_ps;
        case TILEDB_DATETIME_FS:
            return DataType::datetime_fs;
        case TILEDB_DATETIME_AS:
            return DataType::datetime_as;
        case TILEDB_TIME_HR:
            return DataType::time_hr;
        case TILEDB_TIME_MIN:
            return DataType::time_min;
        case TILEDB_TIME_SEC:
            return DataType::time_sec;
        case TILEDB_TIME_MS:
            return DataType::time_ms;
        case TILEDB_TIME_US:
            return DataType::time_us;
        case TILEDB_TIME_NS:
            return DataType::time_ns;
        case TILEDB_TIME_PS:
            return DataType::time_ps;
        case TILEDB_TIME_FS:
            return DataType::time_fs;
        case TILEDB_TIME_AS:
            return DataType::time_as;
        case TILEDB_CHAR:
            return DataType::character;
        case TILEDB_STRING_ASCII:
            return DataType::string_ascii;
        case TILEDB_STRING_UTF8:
            return DataType::string_utf8;
        case TILEDB_STRING_UTF16:
            return DataType::string_utf16;
        case TILEDB_STRING_UTF32:
            return DataType::string_utf32;
        case TILEDB_STRING_UCS2:
            return DataType::string_ucs2;
        case TILEDB_STRING_UCS4:
            return DataType::string_ucs4;
        case TILEDB_BLOB:
            return DataType::blob;
        case TILEDB_GEOM_WKB:
            return DataType::geom_wkb;
        case TILEDB_GEOM_WKT:
            return DataType::geom_wkt;
        default:
            throw std::invalid_argument("Unsupported datatype");
    }
}

// DataType arrow_to_datatype(std::string_view type, std::string_view type_metadata);

constexpr tiledb_datatype_t datatype_to_tiledb(DataType type) {
    switch (type) {
        case DataType::boolean:
            return TILEDB_BOOL;
        case DataType::int8:
            return TILEDB_INT8;
        case DataType::int16:
            return TILEDB_INT16;
        case DataType::int32:
            return TILEDB_INT32;
        case DataType::int64:
            return TILEDB_INT64;
        case DataType::uint8:
            return TILEDB_UINT8;
        case DataType::uint16:
            return TILEDB_UINT16;
        case DataType::uint32:
            return TILEDB_UINT32;
        case DataType::uint64:
            return TILEDB_UINT64;
        case DataType::float32:
            return TILEDB_FLOAT32;
        case DataType::float64:
            return TILEDB_FLOAT64;
        case DataType::datetime_year:
            return TILEDB_DATETIME_YEAR;
        case DataType::datetime_month:
            return TILEDB_DATETIME_MONTH;
        case DataType::datetime_week:
            return TILEDB_DATETIME_WEEK;
        case DataType::datetime_day:
            return TILEDB_DATETIME_DAY;
        case DataType::datetime_hr:
            return TILEDB_DATETIME_HR;
        case DataType::datetime_min:
            return TILEDB_DATETIME_MIN;
        case DataType::datetime_sec:
            return TILEDB_DATETIME_SEC;
        case DataType::datetime_ms:
            return TILEDB_DATETIME_MS;
        case DataType::datetime_us:
            return TILEDB_DATETIME_US;
        case DataType::datetime_ns:
            return TILEDB_DATETIME_NS;
        case DataType::datetime_ps:
            return TILEDB_DATETIME_PS;
        case DataType::datetime_fs:
            return TILEDB_DATETIME_FS;
        case DataType::datetime_as:
            return TILEDB_DATETIME_AS;
        case DataType::time_hr:
            return TILEDB_TIME_HR;
        case DataType::time_min:
            return TILEDB_TIME_MIN;
        case DataType::time_sec:
            return TILEDB_TIME_SEC;
        case DataType::time_ms:
            return TILEDB_TIME_MS;
        case DataType::time_us:
            return TILEDB_TIME_US;
        case DataType::time_ns:
            return TILEDB_TIME_NS;
        case DataType::time_ps:
            return TILEDB_TIME_PS;
        case DataType::time_fs:
            return TILEDB_TIME_FS;
        case DataType::time_as:
            return TILEDB_TIME_AS;
        case DataType::character:
            return TILEDB_CHAR;
        case DataType::string_ascii:
            return TILEDB_STRING_ASCII;
        case DataType::string_utf8:
            return TILEDB_STRING_UTF8;
        case DataType::string_utf16:
            return TILEDB_STRING_UTF16;
        case DataType::string_utf32:
            return TILEDB_STRING_UTF32;
        case DataType::string_ucs2:
            return TILEDB_STRING_UCS2;
        case DataType::string_ucs4:
            return TILEDB_STRING_UCS4;
        case DataType::blob:
            return TILEDB_BLOB;
        case DataType::geom_wkb:
            return TILEDB_GEOM_WKB;
        case DataType::geom_wkt:
            return TILEDB_GEOM_WKT;
        default:
            throw std::invalid_argument("Unsupported datatype");
    }
}

constexpr std::string_view datatype_to_arrow(DataType type, bool use_large = true) {
    auto u = use_large ? "U"sv : "u"sv;
    auto z = use_large ? "Z"sv : "z"sv;

    switch (type) {
        case DataType::boolean:
            return "b"sv;
        case DataType::int8:
            return "c"sv;
        case DataType::int16:
            return "s"sv;
        case DataType::int32:
            return "i"sv;
        case DataType::int64:
            return "l"sv;
        case DataType::uint8:
            return "C"sv;
        case DataType::uint16:
            return "S"sv;
        case DataType::uint32:
            return "I"sv;
        case DataType::uint64:
            return "L"sv;
        case DataType::float32:
            return "f"sv;
        case DataType::float64:
            return "g"sv;
        case DataType::datetime_sec:
            return "tss:"sv;
        case DataType::datetime_ms:
            return "tsm:"sv;
        case DataType::datetime_us:
            return "tsu:"sv;
        case DataType::datetime_ns:
            return "tsn:"sv;
        case DataType::character:
            return z;
        case DataType::string_ascii:
            return u;
        case DataType::string_utf8:
            return u;
        case DataType::blob:
            return z;
        case DataType::geom_wkb:
            return z;
        case DataType::geom_wkt:
            return u;
        default:
            throw std::invalid_argument("Unsuported type by Arrow");
    }
}

constexpr std::optional<std::string_view> datatype_to_arrow_metadata(DataType type) {
    switch (type) {
        case DataType::geom_wkb:
            return "WKB"sv;
        case DataType::geom_wkt:
            return "WKT"sv;
        default:
            return std::nullopt;
    }
}
}  // namespace impl

enum class DataTypeFormat { SOMA, TILEDB, ARROW, ARROW_METADATA };

template <DataTypeFormat Format, bool use_large = true>
constexpr auto as(DataType type) {
    if constexpr (Format == DataTypeFormat::SOMA) {
        return type;
    } else if constexpr (Format == DataTypeFormat::TILEDB) {
        return impl::datatype_to_tiledb(type);
    } else if constexpr (Format == DataTypeFormat::ARROW) {
        return impl::datatype_to_arrow(type, use_large);
    } else {
        return impl::datatype_to_arrow_metadata(type);
    }
}

template <DataTypeFormat Format, bool use_large = true>
constexpr auto as(tiledb_datatype_t type) {
    if constexpr (Format == DataTypeFormat::SOMA) {
        return impl::tiledb_to_datatype(type);
    } else if constexpr (Format == DataTypeFormat::TILEDB) {
        return type;
    } else if constexpr (Format == DataTypeFormat::ARROW) {
        return impl::datatype_to_arrow(impl::tiledb_to_datatype(type), use_large);
    } else {
        return impl::datatype_to_arrow_metadata(impl::tiledb_to_datatype(type));
    }
}
}  // namespace tiledbsoma::common::type

#endif