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

SOMAQueryCondition::SOMAQueryCondition(const QueryCondition& qc)
    : qc_{qc} {
}

SOMAQueryCondition::SOMAQueryCondition(QueryCondition&& qc)
    : qc_{qc} {
}

SOMAQueryCondition& SOMAQueryCondition::combine_with_and(const SOMAQueryCondition& other) {
    if (not other.is_initialized()) {
        // No-op.
        return *this;
    }
    if (qc_.has_value()) {
        qc_ = qc_->combine(other.query_condition(), TILEDB_AND);
    } else {
        qc_ = other.query_condition();
    }
    return *this;
}

SOMAQueryCondition& SOMAQueryCondition::combine_with_or(const SOMAQueryCondition& other) {
    if (not other.is_initialized()) {
        // No-op.
        return *this;
    }
    if (qc_.has_value()) {
        qc_ = qc_->combine(other.query_condition(), TILEDB_OR);
    } else {
        qc_ = other.query_condition();
    }
    return *this;
}

SOMACoordQueryCondition::SOMACoordQueryCondition(const SOMAContext& ctx)
    : ctx_{ctx.tiledb_ctx()}
    , qc_{*ctx_}
    , initialized_{false} {
}

}  // namespace tiledbsoma
