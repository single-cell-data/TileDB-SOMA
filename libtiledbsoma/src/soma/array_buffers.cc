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
#include "column_buffer.h"
#include "common/logging/impl/logger.h"

namespace tiledbsoma {

using namespace tiledb;

bool ArrayBuffers::use_memory_pool(const std::shared_ptr<tiledb::Array>& array) {
    bool use_memory_pool = false;
    auto config = array->config();
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

    MemoryMode mode = ColumnBuffer::memory_mode(array.config());
    const tiledb::ArraySchema schema = array.schema();
    const tiledb::Context& context = array.context();
    // Split memory budget to each column depending on the byte size of each columns element
    // Var sized columns will be allocated the same as an 8 byte datatype

    for (const auto& name : names_) {
        if (schema.has_attribute(name)) {
            tiledb::Attribute attribute = schema.attribute(name);

            auto [column_budget, num_cells] = strategy_->get_buffer_sizes(attribute);
            auto enum_name = AttributeExperimental::get_enumeration_name(context, attribute);
            std::optional<Enumeration> enumeration = std::nullopt;
            bool is_ordered = false;
            if (enum_name.has_value()) {
                enumeration = std::make_optional<Enumeration>(
                    ArrayExperimental::get_enumeration(context, array, *enum_name));
                is_ordered = enumeration->ordered();
            }

            buffers_.insert(
                std::make_pair(
                    name,
                    std::make_shared<CArrayColumnBuffer>(
                        name,
                        attribute.type(),
                        num_cells,
                        column_budget,
                        attribute.variable_sized(),
                        attribute.nullable(),
                        enumeration,
                        is_ordered,
                        mode)));
        }
        // Else check if column is a TileDB dimension
        else if (schema.domain().has_dimension(name)) {
            tiledb::Dimension dimension = schema.domain().dimension(name);

            auto [column_budget, num_cells] = strategy_->get_buffer_sizes(dimension);
            bool is_var = dimension.cell_val_num() == TILEDB_VAR_NUM || dimension.type() == TILEDB_STRING_ASCII ||
                          dimension.type() == TILEDB_STRING_UTF8;

            buffers_.insert(
                std::make_pair(
                    name,
                    std::make_shared<CArrayColumnBuffer>(
                        name, dimension.type(), num_cells, column_budget, is_var, false, std::nullopt, false, mode)));
        }
    }
}

void ArrayBuffers::emplace(const std::string& name, std::shared_ptr<ColumnBuffer> buffer) {
    if (contains(name)) {
        throw TileDBSOMAError(fmt::format("[ArrayBuffers] column '{}' already exists", name));
    }
    names_.push_back(name);
    buffers_.emplace(name, buffer);
}

void ArrayBuffers::expand_buffers() {
    for (const auto& [name, buffer] : buffers_) {
        buffer->resize(
            buffer->max_size() * DEFAULT_BUFFER_EXPANSION_FACTOR,
            buffer->max_num_cells() * DEFAULT_BUFFER_EXPANSION_FACTOR);
    }
}
}  // namespace tiledbsoma
