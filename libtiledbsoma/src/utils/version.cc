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
#include "logger.h"
namespace tiledbsoma::version {

std::string as_string() {
    int major, minor, patch;
    tiledb_version(&major, &minor, &patch);
    return std::format("libtiledb={}.{}.{}", major, minor, patch);
}

std::tuple<int, int, int> embedded_version_triple() {
    int major, minor, patch;
    tiledb_version(&major, &minor, &patch);
    return std::make_tuple(major, minor, patch);
}

};  // namespace tiledbsoma::version
