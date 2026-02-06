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

#include <concepts>
#include <numeric>
#include <tiledb/tiledb>

#include "column_buffer.h"
#include "column_buffer_strategies.h"

#include "../logging/impl/logger.h"

namespace tiledbsoma::common {

bool ArrayBuffers::use_memory_pool(const tiledb::Array& array) {
    bool use_memory_pool = false;
    const tiledb::Config config = array.config();
    if (config.contains(CONFIG_KEY_USE_MEMORY_POOL)) {
        use_memory_pool = config.get(CONFIG_KEY_USE_MEMORY_POOL) == "true";
    }

    return use_memory_pool;
}

ArrayBuffers::ArrayBuffers(
    const std::vector<std::string>& names,
    const tiledb::Array& array,
    std::unique_ptr<ColumnBufferAllocationStrategy> strategy)
    : names_(names)
    , strategy_(std::move(strategy)) {
    if (!strategy_) {
        strategy_ = std::make_unique<BasicAllocationStrategy>(array);
    }

    common::MemoryMode mode = common::ColumnBuffer::memory_mode(array.config());

    for (const auto& name : names_) {
        buffers_.insert(std::make_pair(name, common::CArrayColumnBuffer::create(array, name, strategy_.get(), mode)));
    }
}

bool ArrayBuffers::contains(const std::string& name) const {
    return buffers_.find(name) != buffers_.end();
}

void ArrayBuffers::emplace(const std::string& name, std::shared_ptr<common::ColumnBuffer> buffer) {
    if (contains(name)) {
        throw std::runtime_error(fmt::format("[ArrayBuffers] column '{}' already exists", name));
    }
    names_.push_back(name);
    buffers_.emplace(name, buffer);
}

const std::vector<std::string>& ArrayBuffers::names() const {
    return names_;
}

uint64_t ArrayBuffers::num_rows() const {
    return buffers_.at(names_.front())->cell_count();
}

void ArrayBuffers::expand_buffers() {
    for (const auto& name : names_) {
        std::shared_ptr<common::ReadColumnBuffer> buffer = at<common::ReadColumnBuffer>(name);
        buffer->resize(
            buffer->data_capacity() * DEFAULT_BUFFER_EXPANSION_FACTOR,
            buffer->cell_capacity() * DEFAULT_BUFFER_EXPANSION_FACTOR);
    }
}
}  // namespace tiledbsoma::common
