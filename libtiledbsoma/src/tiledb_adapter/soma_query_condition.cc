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
#include "../utils/logger.h"

namespace tiledbsoma {
using namespace tiledb;

/**************************************
 * SOMAQueryCondition
**************************************/

SOMAQueryCondition::SOMAQueryCondition(const QueryCondition& qc)
    : qc_{qc} {
}

SOMAQueryCondition::SOMAQueryCondition(QueryCondition&& qc)
    : qc_{qc} {
}

SOMACoordQueryCondition::SOMACoordQueryCondition(const SOMAContext& ctx, const std::vector<std::string>& dim_names)
    : ctx_{ctx.tiledb_ctx()}
    , dim_names_{dim_names}
    , coord_qc_(dim_names.size()) {
}

bool SOMACoordQueryCondition::is_initialized() const {
    return std::any_of(coord_qc_.cbegin(), coord_qc_.cend(), [](auto qc) { return qc.is_initialized(); });
}

SOMAQueryCondition SOMAQueryCondition::create_from_range(
    const Context& ctx, const std::string& column_name, const std::string& start_value, const std::string& stop_value) {
    if (stop_value < start_value) {
        throw std::invalid_argument(
            fmt::format(
                "Cannot set range [{}, {}] on column '{}'. Invalid range: the final value must be greater "
                "than or equal to the starting value.",
                start_value,
                stop_value,
                column_name));
    }
    return SOMAQueryCondition(
        QueryCondition::create(ctx, column_name, start_value, TILEDB_GE)
            .combine(QueryCondition::create(ctx, column_name, stop_value, TILEDB_LE), TILEDB_AND));
}

SOMAQueryCondition SOMAQueryCondition::create_from_points(
    const Context& ctx, const std::string& column_name, std::span<std::string> values) {
    if (values.empty()) {
        throw std::invalid_argument(
            fmt::format("Cannot set coordinates on column '{}'. No coordinates provided.", column_name));
    }
    // Using C API because C++ API only supports std::vector, not std::span.
    uint64_t data_size = 0;
    for (auto& val : values) {
        data_size += val.size();
    }
    std::vector<uint8_t> data(data_size);
    std::vector<uint64_t> offsets{};
    uint64_t curr_offset = 0;
    for (auto& val : values) {
        offsets.push_back(curr_offset);
        memcpy(data.data() + curr_offset, val.data(), val.size());
        curr_offset += val.size();
    }
    tiledb_query_condition_t* qc;
    ctx.handle_error(tiledb_query_condition_alloc_set_membership(
        ctx.ptr().get(),
        column_name.c_str(),
        data.data(),
        data.size(),
        offsets.data(),
        offsets.size() * sizeof(uint64_t),
        TILEDB_IN,
        &qc));
    return QueryCondition(ctx, qc);
}

/**************************************
 * SOMACoordQueryCondition
**************************************/

SOMACoordQueryCondition& SOMACoordQueryCondition::add_column_selection(
    int64_t col_index, SOMAColumnSelection<std::string> selection) {
    std::visit(
        [&](auto&& val) {
            using S = std::decay_t<decltype(val)>;
            if constexpr (std::is_same_v<S, SOMASliceSelection<std::string>>) {
                add_coordinate_query_condition(
                    col_index,
                    SOMAQueryCondition::create_from_range(*ctx_, dim_names_[col_index], val.start, val.stop));

            } else if constexpr (std::is_same_v<S, SOMAPointSelection<std::string>>) {
                add_coordinate_query_condition(
                    col_index, SOMAQueryCondition::create_from_points(*ctx_, dim_names_[col_index], val.points));
            }
            // Otherwise monostate: do nothing.
        },
        selection);
    return *this;
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
