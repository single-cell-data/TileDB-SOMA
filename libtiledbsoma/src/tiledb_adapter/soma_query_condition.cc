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
 * SOMAValueFilter
**************************************/

SOMAValueFilter::SOMAValueFilter(const QueryCondition& qc)
    : qc_{qc} {
}

SOMAValueFilter::SOMAValueFilter(QueryCondition&& qc)
    : qc_{qc} {
}

SOMAValueFilter SOMAValueFilter::create_from_slice(
    const Context& ctx, const std::string& column_name, const SOMASliceSelection<std::string>& slice) {
    if (slice.start.has_value() && slice.stop.has_value()) {
        return SOMAValueFilter(
            QueryCondition::create(ctx, column_name, slice.start.value(), TILEDB_GE)
                .combine(QueryCondition::create(ctx, column_name, slice.stop.value(), TILEDB_LE), TILEDB_AND));
    }
    if (slice.start.has_value()) {
        return SOMAValueFilter(QueryCondition::create(ctx, column_name, slice.start.value(), TILEDB_GE));
    }
    if (slice.stop.has_value()) {
        return SOMAValueFilter(QueryCondition::create(ctx, column_name, slice.stop.value(), TILEDB_LE));
    }
    return SOMAValueFilter();
}

SOMAValueFilter SOMAValueFilter::create_from_points(
    const Context& ctx, const std::string& column_name, SOMAPointSelection<std::string> values) {
    if (values.points.empty()) {
        throw std::invalid_argument(
            fmt::format("Cannot set coordinates on column '{}'. No coordinates provided.", column_name));
    }
    // Using C API because C++ API only supports std::vector, not std::span.
    uint64_t data_size = 0;
    for (auto& val : values.points) {
        data_size += val.size();
    }
    std::vector<uint8_t> data(data_size);
    std::vector<uint64_t> offsets{};
    uint64_t curr_offset = 0;
    for (auto& val : values.points) {
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
 * CoordinateValueFilter
**************************************/

CoordinateValueFilter::CoordinateValueFilter(
    std::shared_ptr<Array> array,
    std::shared_ptr<SOMAContext> ctx,
    std::vector<std::shared_ptr<SOMAColumn>> index_columns,
    Domainish domain_kind)
    : ctx_{ctx}
    , array_{array}
    , index_columns_{index_columns}
    , domain_kind_{domain_kind}
    , coord_qc_(index_columns_.size()) {
}

bool CoordinateValueFilter::is_initialized() const {
    return std::any_of(coord_qc_.cbegin(), coord_qc_.cend(), [](auto qc) { return qc.is_initialized(); });
}

CoordinateValueFilter& CoordinateValueFilter::add_column_selection(
    int64_t col_index, SOMAColumnSelection<std::string> selection) {
    std::visit(
        [&](auto&& val) {
            using S = std::decay_t<decltype(val)>;
            if constexpr (std::is_same_v<S, SOMASliceSelection<std::string>>) {
                add_slice(col_index, val);
            } else if constexpr (std::is_same_v<S, SOMAPointSelection<std::string>>) {
                add_points(col_index, val);
            }
            // Otherwise monostate: do nothing.
        },
        selection);
    return *this;
}

SOMAValueFilter CoordinateValueFilter::get_value_filter() const {
    return std::reduce(coord_qc_.cbegin(), coord_qc_.cend(), SOMAValueFilter(), [](const auto& qc1, const auto& qc2) {
        if (qc1.is_initialized()) {
            if (qc2.is_initialized()) {
                return SOMAValueFilter(qc1.query_condition().combine(qc2.query_condition(), TILEDB_AND));
            }
            return qc1;
        }
        return qc2;
    });
}

CoordinateValueFilter& CoordinateValueFilter::add_coordinate_query_condition(int64_t index, SOMAValueFilter&& qc) {
    if (!qc.is_initialized()) {
        // No-op.
        return *this;
    }
    auto& current = coord_qc_[index];
    if (current.is_initialized()) {
        coord_qc_[index] = SOMAValueFilter(current.query_condition().combine(qc.query_condition(), TILEDB_OR));
    } else {
        coord_qc_[index] = qc;
    }
    return *this;
}

}  // namespace tiledbsoma
