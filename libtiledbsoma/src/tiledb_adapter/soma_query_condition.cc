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

/**Throws an error if using a string on a non-string column. */
void CoordinateValueFilter::validate_string_column(std::shared_ptr<SOMAColumn> column) const {
    if (column->domain_type().value() != TILEDB_STRING_ASCII) {
        std::stringstream ss;
        auto tiledb_type = tiledb::impl::type_to_str(column->domain_type().value());
        ss << "Invalid coordinate on column '" << column->name() << "'. Cannot set string value on column with type "
           << tiledb_type << ".";
        throw std::invalid_argument(ss.str());
    }
}

}  // namespace tiledbsoma
