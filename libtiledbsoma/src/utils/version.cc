/**
 * @file   version.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This exposes the version of the TileDB Embedded library in use.
 */

#include "version.h"
#include <tiledb/tiledb>

#include "common/logging/impl/logger.h"

namespace tiledbsoma::version {

std::string as_string() {
    int major, minor, patch;
    tiledb_version(&major, &minor, &patch);
    return fmt::format("libtiledb={}.{}.{}", major, minor, patch);
}

std::tuple<int, int, int> embedded_version_triple() {
    int major, minor, patch;
    tiledb_version(&major, &minor, &patch);
    return std::make_tuple(major, minor, patch);
}

/**
 * @brief Link-time version of TileDB. Used for error checking at runtime.
 *
 * @return std::tuple<int, int, int>
 */
std::tuple<int, int, int> expected_version() {
    return std::make_tuple(TILEDB_VERSION_MAJOR, TILEDB_VERSION_MINOR, TILEDB_VERSION_PATCH);
}

};  // namespace tiledbsoma::version
