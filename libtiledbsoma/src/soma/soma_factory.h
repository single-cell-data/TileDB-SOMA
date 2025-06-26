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

#include <tiledb/tiledb>

#include "../utils/common.h"

using namespace tiledb;
namespace tiledbsoma {

class SOMAContext;

std::string get_soma_type(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp = std::nullopt);

std::string get_tiledb_array_soma_type(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp = std::nullopt);

std::string get_tiledb_group_soma_type(
    std::string_view uri, const SOMAContext& ctx, std::optional<TimestampRange> timestamp = std::nullopt);

}  // namespace tiledbsoma

#endif
