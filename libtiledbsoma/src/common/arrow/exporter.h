/**
 * @file   exporter.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 */

#ifndef COMMON_ARROW_EXPORTER_H
#define COMMON_ARROW_EXPORTER_H

#include <future>
#include <span>
#include <unordered_map>

#include "../query/column_buffer.h"
#include "utils.h"

namespace tiledbsoma::common::arrow {

ArrowTable column_buffer_to_arrow(
    ColumnBuffer* column_buffer,
    const std::unordered_map<std::string, std::shared_future<std::shared_ptr<ArrowBuffer>>>& enumerations,
    bool downcast_dict_of_large_var);

}  // namespace tiledbsoma::common::arrow

#endif