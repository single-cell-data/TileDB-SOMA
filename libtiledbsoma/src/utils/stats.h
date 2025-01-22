/**
 * @file   stats.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This declares the common functions in the API
 */

#ifndef TILEDBSOMA_STATS_H
#define TILEDBSOMA_STATS_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <string>

namespace tiledbsoma::stats {

void enable();
void disable();
void reset();
std::string dump();

};  // namespace tiledbsoma::stats

#endif  // TILEDBSOMA_STATS_H
