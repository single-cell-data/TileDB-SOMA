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
#include "column_buffer.h"

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

ArrayBuffers::ArrayBuffers(const std::vector<std::string>& names, const tiledb::Array& array) {
    size_t memory_budget = DEFAULT_ALLOC_BYTES;
    auto config = array.config();
    if (config.contains(CONFIG_KEY_MEMORY_BUDGET)) {
        auto value_str = config.get(CONFIG_KEY_MEMORY_BUDGET);
        try {
            memory_budget = std::stoull(value_str);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(
                fmt::format(
                    "[ArrayBuffers] Error parsing {}: '{}' ({})", CONFIG_KEY_MEMORY_BUDGET, value_str, e.what()));
        }
    }

    size_t factor = 1;
    if (config.contains(CONFIG_KEY_VAR_SIZED_FACTOR)) {
        auto value_str = config.get(CONFIG_KEY_VAR_SIZED_FACTOR);
        try {
            factor = std::stoull(value_str);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(
                fmt::format(
                    "[ArrayBuffers] Error parsing {}: '{}' ({})", CONFIG_KEY_VAR_SIZED_FACTOR, value_str, e.what()));
        }
    }

    ArraySchema schema = array.schema();
    // Split memory budget to each column depending on the byte size of each columns element
    // Var sized columns will be allocated the same as an 8 byte datatype

    // Ensure minimum buffer size is multiple of 8
    size_t memory_budget_unit = (memory_budget /
                                 std::transform_reduce(
                                     names.begin(),
                                     names.end(),
                                     0L,
                                     std::plus{},
                                     [&](auto name) {
                                         size_t weight = 0;

                                         // Check if column is a TileDB attribute
                                         if (schema.has_attribute(name)) {
                                             Attribute attr = schema.attribute(name);

                                             if (!attr.variable_sized() && attr.cell_val_num() != 1) {
                                                 throw TileDBSOMAError(
                                                     "[ArrayBuffers] Values per cell > 1 is not supported: " + name);
                                             }

                                             weight += attr.nullable() ? 1 : 0;
                                             // If column has variable size add the offset array in the column budget
                                             weight += attr.variable_sized() ? sizeof(uint64_t) * (1 + 2 * factor) :
                                                                               tiledb::impl::type_size(attr.type());

                                             return weight;
                                         }
                                         // Else check if column is a TileDB dimension
                                         else if (schema.domain().has_dimension(name)) {
                                             Dimension dim = schema.domain().dimension(name);

                                             bool is_var = dim.cell_val_num() == TILEDB_VAR_NUM ||
                                                           dim.type() == TILEDB_STRING_ASCII ||
                                                           dim.type() == TILEDB_STRING_UTF8;

                                             if (!is_var && dim.cell_val_num() != 1) {
                                                 throw TileDBSOMAError(
                                                     "[ArrayBuffers] Values per cell > 1 is not supported: " + name);
                                             }

                                             weight += (dim.type() == TILEDB_STRING_ASCII ||
                                                        dim.type() == TILEDB_STRING_UTF8) ?
                                                           sizeof(uint64_t) * (1 + 2 * factor) :
                                                           tiledb::impl::type_size(dim.type());

                                             return weight;
                                         }

                                         throw TileDBSOMAError(
                                             fmt::format("[ArrayBuffers] Missing column name '{}'", name));
                                     }) /
                                 8) *
                                8;

    for (const auto& name : names) {
        names_.push_back(name);

        if (schema.has_attribute(name)) {
            Attribute attr = schema.attribute(name);

            size_t column_budget = (attr.variable_sized() ? sizeof(uint64_t) * 2 * factor :
                                                            tiledb::impl::type_size(attr.type())) *
                                   memory_budget_unit;

            size_t num_cells = memory_budget_unit;

            auto enum_name = AttributeExperimental::get_enumeration_name(schema.context(), attr);
            std::optional<Enumeration> enumeration = std::nullopt;
            bool is_ordered = false;
            if (enum_name.has_value()) {
                auto enmr = ArrayExperimental::get_enumeration(schema.context(), array, *enum_name);
                is_ordered = enmr.ordered();
                enumeration = std::make_optional<Enumeration>(enmr);
            }

            buffers_.insert(
                std::make_pair(
                    name,
                    std::make_shared<CArrayColumnBuffer>(
                        name,
                        attr.type(),
                        num_cells,
                        column_budget,
                        attr.variable_sized(),
                        attr.nullable(),
                        enumeration,
                        is_ordered)));
        }
        // Else check if column is a TileDB dimension
        else if (schema.domain().has_dimension(name)) {
            Dimension dim = schema.domain().dimension(name);

            bool is_var = dim.cell_val_num() == TILEDB_VAR_NUM || dim.type() == TILEDB_STRING_ASCII ||
                          dim.type() == TILEDB_STRING_UTF8;

            // Ensure buffer size is multiple of 8
            size_t column_budget = (is_var ? sizeof(uint64_t) * 2 * factor : tiledb::impl::type_size(dim.type())) *
                                   memory_budget_unit;
            size_t num_cells = memory_budget_unit;

            buffers_.insert(
                std::make_pair(
                    name,
                    std::make_shared<CArrayColumnBuffer>(
                        name, dim.type(), num_cells, column_budget, is_var, false, std::nullopt, false)));
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

}  // namespace tiledbsoma
