/**
 * @file   version.h
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

#ifndef TILEDBSOMA_VERSION_H
#define TILEDBSOMA_VERSION_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <string>

namespace tiledbsoma::version {

std::string as_string();
std::tuple<int, int, int> embedded_version_triple();

};  // namespace tiledbsoma::version

#endif  // TILEDBSOMA_VERSION_H
