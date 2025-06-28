/**
 * @file   soma_factory.h
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

#ifndef TILEDBSOMA_FACTORY_H
#define TILEDBSOMA_FACTORY_H

#include <optional>
#include <string>
#include "../utils/common.h"

using namespace tiledb;
namespace tiledbsoma {

class SOMAContext;

/**
 * Read SOMA datatype metadata from a TileDB object.
 *
 * Warning: This does not validate a SOMA object or metadata beyond the
 * existance and stored type of the SOMA datatype metadata.
 *
 * @param uri The URI for the TileDB array or group.
 * @param ctx A SOMA context object to use for accesing the array or group.
 * @param timestamp An optional timestamp value for time travel.
 */
std::string get_soma_type_metadata_value(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp = std::nullopt);

/**
 * Read SOMA datatype metadata from a TileDB array.
 *
 * Warning: This does not validate a SOMA object or metadata beyond the
 * existance and stored type of the SOMA datatype metadata.
 *
 * @param uri The URI for the TileDB array.
 * @param ctx A SOMA context object to use for accesing the array.
 * @param timestamp An optional timestamp value for time travel.
 */
std::string get_soma_type_metadata_value_from_array(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp = std::nullopt);

/**
 * Read SOMA datatype metadata from a TileDB group.
 *
 * Warning: This does not validate a SOMA object or metadata beyond the
 * existance and stored type of the SOMA datatype metadata.
 *
 * @param uri The URI for the TileDB group.
 * @param ctx A SOMA context object to use for accesing the group.
 * @param timestamp An optional timestamp value for time travel.
 */
std::string get_soma_type_metadata_value_from_group(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp = std::nullopt);

}  // namespace tiledbsoma

#endif
