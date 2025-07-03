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

#include <numeric>
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

SOMACoordQueryCondition::SOMACoordQueryCondition(const SOMAContext& ctx, const std::vector<std::string>& dim_names)
    : ctx_{ctx.tiledb_ctx()}
    , dim_names_{dim_names}
    , coord_qc_(dim_names.size()) {
}

bool SOMACoordQueryCondition::is_initialized() const {
    return std::any_of(coord_qc_.cbegin(), coord_qc_.cend(), [](auto qc) { return qc.is_initialized(); });
}

SOMAQueryCondition SOMACoordQueryCondition::get_soma_query_condition() const {
    return std::reduce(
        coord_qc_.cbegin(), coord_qc_.cend(), SOMAQueryCondition(), [](const auto& qc1, const auto& qc2) {
            if (qc1.is_initialized()) {
                if (qc2.is_initialized()) {
                    return SOMAQueryCondition(qc1.query_condition().combine(qc2.query_condition(), TILEDB_AND));
                }
                return qc1;
            }
            return qc2;
        });
    return SOMAQueryCondition();
}

SOMACoordQueryCondition& SOMACoordQueryCondition::add_coordinate_query_condition(
    int64_t index, SOMAQueryCondition&& qc) {
    if (!qc.is_initialized()) {
        // No-op.
        return *this;
    }
    auto& current = coord_qc_[index];
    if (current.is_initialized()) {
        coord_qc_[index] = SOMAQueryCondition(current.query_condition().combine(qc.query_condition(), TILEDB_OR));
    } else {
        coord_qc_[index] = qc;
    }
    return *this;
}

}  // namespace tiledbsoma
