/**
 * @file   array_buffers.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the ArrayBuffers class.
 */

#include "array_buffers.h"
#include "../utils/logger.h"

#include <format>

namespace tiledbsoma {

using namespace tiledb;

std::shared_ptr<ColumnBuffer> ArrayBuffers::at(const std::string& name) {
    if (!contains(name)) {
        throw TileDBSOMAError(
            std::format("[ArrayBuffers] column '{}' does not exist", name));
    }
    return buffers_[name];
}

void ArrayBuffers::emplace(
    const std::string& name, std::shared_ptr<ColumnBuffer> buffer) {
    if (contains(name)) {
        throw TileDBSOMAError(
            std::format("[ArrayBuffers] column '{}' already exists", name));
    }
    names_.push_back(name);
    buffers_.emplace(name, buffer);
}

}  // namespace tiledbsoma
