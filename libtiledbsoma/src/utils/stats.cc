/**
 * @file   stats.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file provides access to stats from libtiledbsoma's dependency on
 * TileDB Embedded.
 */

#include "utils/stats.h"
#include <tiledb/tiledb>

namespace tiledbsoma::stats {

void enable() {
    tiledb::Stats::enable();
}

void disable() {
    tiledb::Stats::disable();
}

void reset() {
    tiledb::Stats::reset();
}

std::string dump() {
    std::string stats;
    tiledb::Stats::raw_dump(&stats);
    return stats;
}

};  // namespace tiledbsoma::stats
