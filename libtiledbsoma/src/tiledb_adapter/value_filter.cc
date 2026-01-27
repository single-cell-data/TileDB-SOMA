/**
 * @file   value_filter.cc
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

#include "value_filter.h"

#include <numeric>

namespace tiledbsoma {
using namespace tiledb;

ValueFilter::ValueFilter(const QueryCondition& qc)
    : qc_{qc} {
}

ValueFilter::ValueFilter(QueryCondition&& qc)
    : qc_{qc} {
}

}  // namespace tiledbsoma
