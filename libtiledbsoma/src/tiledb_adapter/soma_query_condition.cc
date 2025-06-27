/**
 * @file   soma_query_condition.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines helper functions and classes for using query conditions in SOMA.
 */

#include "soma_query_condition.h"

#include "../utils/common.h"

namespace tiledbsoma {
using namespace tiledb;

SOMACoordQueryCondition::SOMACoordQueryCondition(const SOMAContext& ctx)
    : ctx_{ctx.tiledb_ctx()}
    , qc_{*ctx_}
    , initialized_{false} {
}

}  // namespace tiledbsoma
