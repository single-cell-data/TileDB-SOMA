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
            return DataType::BOOL;
        case TILEDB_INT8:
            return DataType::INT8;
        case TILEDB_INT16:
            return DataType::INT16;
        case TILEDB_INT32:
            return DataType::INT32;
        case TILEDB_INT64:
            return DataType::INT64;
        case TILEDB_UINT8:
            return DataType::UINT8;
        case TILEDB_UINT16:
            return DataType::UINT16;
        case TILEDB_UINT32:
            return DataType::UINT32;
        case TILEDB_UINT64:
            return DataType::UINT64;
        case TILEDB_FLOAT32:
            return DataType::FLOAT32;
        case TILEDB_FLOAT64:
            return DataType::FLOAT64;
        case TILEDB_DATETIME_YEAR:
            return DataType::DATETIME_YEAR;
        case TILEDB_DATETIME_MONTH:
            return DataType::DATETIME_MONTH;
        case TILEDB_DATETIME_WEEK:
            return DataType::DATETIME_WEEK;
        case TILEDB_DATETIME_DAY:
            return DataType::DATETIME_DAY;
        case TILEDB_DATETIME_HR:
            return DataType::DATETIME_HR;
        case TILEDB_DATETIME_MIN:
            return DataType::DATETIME_MIN;
        case TILEDB_DATETIME_SEC:
            return DataType::DATETIME_SEC;
        case TILEDB_DATETIME_MS:
            return DataType::DATETIME_MS;
        case TILEDB_DATETIME_US:
            return DataType::DATETIME_US;
        case TILEDB_DATETIME_NS:
            return DataType::DATETIME_NS;
        case TILEDB_DATETIME_PS:
            return DataType::DATETIME_PS;
        case TILEDB_DATETIME_FS:
            return DataType::DATETIME_FS;
        case TILEDB_DATETIME_AS:
            return DataType::DATETIME_AS;
        case TILEDB_TIME_HR:
            return DataType::TIME_HR;
        case TILEDB_TIME_MIN:
            return DataType::TIME_MIN;
        case TILEDB_TIME_SEC:
            return DataType::TIME_SEC;
        case TILEDB_TIME_MS:
            return DataType::TIME_MS;
        case TILEDB_TIME_US:
            return DataType::TIME_US;
        case TILEDB_TIME_NS:
            return DataType::TIME_NS;
        case TILEDB_TIME_PS:
            return DataType::TIME_PS;
        case TILEDB_TIME_FS:
            return DataType::TIME_FS;
        case TILEDB_TIME_AS:
            return DataType::TIME_AS;
        case TILEDB_CHAR:
            return DataType::CHAR;
        case TILEDB_STRING_ASCII:
            return DataType::STRING_ASCII;
        case TILEDB_STRING_UTF8:
            return DataType::STRING_UTF8;
        case TILEDB_STRING_UTF16:
            return DataType::STRING_UTF16;
        case TILEDB_STRING_UTF32:
            return DataType::STRING_UTF32;
        case TILEDB_STRING_UCS2:
            return DataType::STRING_UCS2;
        case TILEDB_STRING_UCS4:
            return DataType::STRING_UCS4;
        case TILEDB_BLOB:
            return DataType::BLOB;
        case TILEDB_GEOM_WKB:
            return DataType::GEOM_WKB;
        case TILEDB_GEOM_WKT:
            return DataType::GEOM_WKT;
        default:
            throw std::invalid_argument("Unsupported datatype");
    }
}

// DataType arrow_to_datatype(std::string_view type, std::string_view type_metadata);

constexpr tiledb_datatype_t datatype_to_tiledb(DataType type) {
    switch (type) {
        case DataType::BOOL:
            return TILEDB_BOOL;
        case DataType::INT8:
            return TILEDB_INT8;
        case DataType::INT16:
            return TILEDB_INT16;
        case DataType::INT32:
            return TILEDB_INT32;
        case DataType::INT64:
            return TILEDB_INT64;
        case DataType::UINT8:
            return TILEDB_UINT8;
        case DataType::UINT16:
            return TILEDB_UINT16;
        case DataType::UINT32:
            return TILEDB_UINT32;
        case DataType::UINT64:
            return TILEDB_UINT64;
        case DataType::FLOAT32:
            return TILEDB_FLOAT32;
        case DataType::FLOAT64:
            return TILEDB_FLOAT64;
        case DataType::DATETIME_YEAR:
            return TILEDB_DATETIME_YEAR;
        case DataType::DATETIME_MONTH:
            return TILEDB_DATETIME_MONTH;
        case DataType::DATETIME_WEEK:
            return TILEDB_DATETIME_WEEK;
        case DataType::DATETIME_DAY:
            return TILEDB_DATETIME_DAY;
        case DataType::DATETIME_HR:
            return TILEDB_DATETIME_HR;
        case DataType::DATETIME_MIN:
            return TILEDB_DATETIME_MIN;
        case DataType::DATETIME_SEC:
            return TILEDB_DATETIME_SEC;
        case DataType::DATETIME_MS:
            return TILEDB_DATETIME_MS;
        case DataType::DATETIME_US:
            return TILEDB_DATETIME_US;
        case DataType::DATETIME_NS:
            return TILEDB_DATETIME_NS;
        case DataType::DATETIME_PS:
            return TILEDB_DATETIME_PS;
        case DataType::DATETIME_FS:
            return TILEDB_DATETIME_FS;
        case DataType::DATETIME_AS:
            return TILEDB_DATETIME_AS;
        case DataType::TIME_HR:
            return TILEDB_TIME_HR;
        case DataType::TIME_MIN:
            return TILEDB_TIME_MIN;
        case DataType::TIME_SEC:
            return TILEDB_TIME_SEC;
        case DataType::TIME_MS:
            return TILEDB_TIME_MS;
        case DataType::TIME_US:
            return TILEDB_TIME_US;
        case DataType::TIME_NS:
            return TILEDB_TIME_NS;
        case DataType::TIME_PS:
            return TILEDB_TIME_PS;
        case DataType::TIME_FS:
            return TILEDB_TIME_FS;
        case DataType::TIME_AS:
            return TILEDB_TIME_AS;
        case DataType::CHAR:
            return TILEDB_CHAR;
        case DataType::STRING_ASCII:
            return TILEDB_STRING_ASCII;
        case DataType::STRING_UTF8:
            return TILEDB_STRING_UTF8;
        case DataType::STRING_UTF16:
            return TILEDB_STRING_UTF16;
        case DataType::STRING_UTF32:
            return TILEDB_STRING_UTF32;
        case DataType::STRING_UCS2:
            return TILEDB_STRING_UCS2;
        case DataType::STRING_UCS4:
            return TILEDB_STRING_UCS4;
        case DataType::BLOB:
            return TILEDB_BLOB;
        case DataType::GEOM_WKB:
            return TILEDB_GEOM_WKB;
        case DataType::GEOM_WKT:
            return TILEDB_GEOM_WKT;
        default:
            throw std::invalid_argument("Unsupported datatype");
    }
}

constexpr std::string_view datatype_to_arrow(DataType type, bool use_large = true) {
    auto u = use_large ? "U"sv : "u"sv;
    auto z = use_large ? "Z"sv : "z"sv;

    switch (type) {
        case DataType::BOOL:
            return "b"sv;
        case DataType::INT8:
            return "c"sv;
        case DataType::INT16:
            return "s"sv;
        case DataType::INT32:
            return "i"sv;
        case DataType::INT64:
            return "l"sv;
        case DataType::UINT8:
            return "C"sv;
        case DataType::UINT16:
            return "S"sv;
        case DataType::UINT32:
            return "I"sv;
        case DataType::UINT64:
            return "L"sv;
        case DataType::FLOAT32:
            return "f"sv;
        case DataType::FLOAT64:
            return "g"sv;
        case DataType::DATETIME_SEC:
            return "tss:"sv;
        case DataType::DATETIME_MS:
            return "tsm:"sv;
        case DataType::DATETIME_US:
            return "tsu:"sv;
        case DataType::DATETIME_NS:
            return "tsn:"sv;
        case DataType::CHAR:
            return z;
        case DataType::STRING_ASCII:
            return u;
        case DataType::STRING_UTF8:
            return u;
        case DataType::BLOB:
            return z;
        case DataType::GEOM_WKB:
            return z;
        case DataType::GEOM_WKT:
            return u;
        default:
            throw std::invalid_argument("Unsuported type by Arrow");
    }
}

constexpr std::optional<std::string_view> datatype_to_arrow_metadata(DataType type) {
    switch (type) {
        case DataType::GEOM_WKB:
            return "WKB"sv;
        case DataType::GEOM_WKT:
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