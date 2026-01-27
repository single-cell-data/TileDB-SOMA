/**
 * @file   soma_factory.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * Methods for accessing SOMAObjects.
 *
 */

#include "soma_factory.h"

#include <tiledb/tiledb>
#include "common/logging/impl/logger.h"
#include "soma_context.h"

namespace tiledbsoma {
using namespace tiledb;

std::string get_soma_type_metadata_value(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp) {
    auto tiledb_type = Object::object(*ctx.tiledb_ctx(), std::string(uri)).type();

    switch (tiledb_type) {
        case Object::Type::Array:
            return get_soma_type_metadata_value_from_array(uri, ctx, timestamp);
        case Object::Type::Group:
            return get_soma_type_metadata_value_from_group(uri, ctx, timestamp);
        default:
            throw TileDBSOMAError(fmt::format("The object at URI '{}' is not a valid TileDB array or group.", uri));
    }
}

std::string get_soma_type_metadata_value_from_array(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp) {
    auto temporal_policy = timestamp.has_value() ?
                               TemporalPolicy(TimestampStartEnd, timestamp->first, timestamp->second) :
                               TemporalPolicy();
    Array array(*ctx.tiledb_ctx(), std::string(uri), TILEDB_READ, temporal_policy);

    // Get and return datatype.
    tiledb_datatype_t value_type{};
    auto has_soma_type = array.has_metadata(SOMA_OBJECT_TYPE_KEY, &value_type);
    if (!has_soma_type) {
        throw TileDBSOMAError(
            fmt::format(
                "Cannot get the SOMA type of the TileDB array at URI '{}'. Missing the required metadata key '{}' for "
                "the SOMA type.",
                uri,
                SOMA_OBJECT_TYPE_KEY));
    }
    if (value_type != TILEDB_STRING_UTF8) {
        const char* expected_type = nullptr;
        tiledb_datatype_to_str(TILEDB_STRING_UTF8, &expected_type);
        const char* actual_type = nullptr;
        tiledb_datatype_to_str(value_type, &actual_type);
        throw TileDBSOMAError(
            fmt::format(
                "Cannot get the SOMA type of the TileDB array at URI '{}'. Metadata for '{}' has datatype '{}'. "
                "Expected datatype '{}'.",
                uri,
                SOMA_OBJECT_TYPE_KEY,
                std::string(actual_type),
                std::string(expected_type)));
    }
    uint32_t value_num{};
    const void* value{};
    array.get_metadata(SOMA_OBJECT_TYPE_KEY, &value_type, &value_num, &value);
    return std::string(static_cast<const char*>(value), value_num);
}

std::string get_soma_type_metadata_value_from_group(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp) {
    auto cfg = ctx.tiledb_ctx()->config();
    if (timestamp.has_value()) {
        if (timestamp->first > timestamp->second) {
            throw std::invalid_argument(
                fmt::format(
                    "Invalid timestamp range. Timestamp start is greater than timestamp end. Timestamp start={}, "
                    "timestamp end={}.",
                    timestamp->first,
                    timestamp->second));
        }
        cfg["sm.group.timestamp_start"] = timestamp->first;
        cfg["sm.group.timestamp_end"] = timestamp->second;
    }
    Group group(*ctx.tiledb_ctx(), std::string(uri), TILEDB_READ, cfg);

    // Get and return datatype.
    tiledb_datatype_t value_type{};
    auto has_soma_type = group.has_metadata(SOMA_OBJECT_TYPE_KEY, &value_type);
    if (!has_soma_type) {
        throw TileDBSOMAError(
            fmt::format(
                "Cannot get the SOMA type of the TileDB group at URI '{}'. Missing the required metadata key '{}' for "
                "the SOMA type.",
                uri,
                SOMA_OBJECT_TYPE_KEY));
    }
    if (value_type != TILEDB_STRING_UTF8) {
        const char* expected_type = nullptr;
        tiledb_datatype_to_str(TILEDB_STRING_UTF8, &expected_type);
        const char* actual_type = nullptr;
        tiledb_datatype_to_str(value_type, &actual_type);
        throw TileDBSOMAError(
            fmt::format(
                "Cannot get the SOMA type of the TileDB group at URI '{}'. Metadata for '{}' has datatype '{}'. "
                "Expected datatype '{}'.",
                uri,
                SOMA_OBJECT_TYPE_KEY,
                std::string(actual_type),
                std::string(expected_type)));
    }
    uint32_t value_num{};
    const void* value{};
    group.get_metadata(SOMA_OBJECT_TYPE_KEY, &value_type, &value_num, &value);
    return std::string(static_cast<const char*>(value), value_num);
}

}  // namespace tiledbsoma
