#include "column_buffer_strategies.h"

#include "../logging/impl/logger.h"

#include <stdexcept>

namespace tiledbsoma {

BasicAllocationStrategy::BasicAllocationStrategy(const tiledb::Array& array) {
    const tiledb::Config config = array.config();

    if (config.contains(CONFIG_KEY_INIT_BYTES)) {
        auto value_str = config.get(CONFIG_KEY_INIT_BYTES.data());
        try {
            buffer_size = std::stoull(value_str);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                fmt::format(
                    "[BasicAllocationStrategy] Error parsing {}: '{}' ({})",
                    CONFIG_KEY_INIT_BYTES,
                    value_str,
                    e.what()));
        }
    }
}

std::pair<size_t, size_t> BasicAllocationStrategy::get_buffer_sizes(
    const std::variant<tiledb::Attribute, tiledb::Dimension> column) const {
    return std::visit(
        [&](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;

            bool is_var = false;
            if constexpr (std::is_same_v<T, tiledb::Attribute>) {
                is_var = arg.variable_sized();
            } else {
                is_var = arg.cell_val_num() == TILEDB_VAR_NUM || arg.type() == TILEDB_STRING_ASCII ||
                         arg.type() == TILEDB_STRING_UTF8;
            }

            size_t num_cells = is_var ? buffer_size / sizeof(uint64_t) :
                                        buffer_size / tiledb::impl::type_size(arg.type());

            return std::make_pair(buffer_size, num_cells);
        },
        column);
}

MemoryPoolAllocationStrategy::MemoryPoolAllocationStrategy(std::span<std::string> columns, const tiledb::Array& array) {
    const tiledb::Config config = array.config();
    const tiledb::ArraySchema schema = array.schema();

    size_t memory_budget = DEFAULT_MEMORY_POOL;
    if (config.contains(CONFIG_KEY_MEMORY_POOL_SIZE)) {
        auto value_str = config.get(CONFIG_KEY_MEMORY_POOL_SIZE.data());
        try {
            memory_budget = std::stoull(value_str);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                fmt::format(
                    "[MemoryPoolAllocationStrategy] Error parsing {}: '{}' ({})",
                    CONFIG_KEY_MEMORY_POOL_SIZE,
                    value_str,
                    e.what()));
        }
    }

    var_size_expansion_factor = DEFAULT_VAR_SIZE_FACTOR;
    if (config.contains(CONFIG_KEY_VAR_SIZED_FACTOR)) {
        auto value_str = config.get(CONFIG_KEY_VAR_SIZED_FACTOR.data());
        try {
            var_size_expansion_factor = std::max(1ull, std::stoull(value_str));
        } catch (const std::exception& e) {
            throw std::runtime_error(
                fmt::format(
                    "[MemoryPoolAllocationStrategy] Error parsing {}: '{}' ({})",
                    CONFIG_KEY_VAR_SIZED_FACTOR,
                    value_str,
                    e.what()));
        }
    }

    size_t weight = 0;
    for (const auto& column : columns) {
        if (schema.has_attribute(column)) {
            tiledb::Attribute attr = schema.attribute(column);

            if (!attr.variable_sized() && attr.cell_val_num() != 1) {
                throw std::runtime_error(
                    "[MemoryPoolAllocationStrategy] Values per cell > 1 is not supported: " + column);
            }

            weight += attr.nullable() ? 1 : 0;
            // If column has variable size add the offset array in the column budget
            weight += attr.variable_sized() ? sizeof(uint64_t) * (1 + var_size_expansion_factor) :
                                              tiledb::impl::type_size(attr.type());
        }
        // Else check if column is a TileDB dimension
        else if (schema.domain().has_dimension(column)) {
            tiledb::Dimension dim = schema.domain().dimension(column);

            bool is_var = dim.cell_val_num() == TILEDB_VAR_NUM || dim.type() == TILEDB_STRING_ASCII ||
                          dim.type() == TILEDB_STRING_UTF8;

            if (!is_var && dim.cell_val_num() != 1) {
                throw std::runtime_error(
                    "[MemoryPoolAllocationStrategy] Values per cell > 1 is not supported: " + column);
            }

            weight += (dim.type() == TILEDB_STRING_ASCII || dim.type() == TILEDB_STRING_UTF8) ?
                          sizeof(uint64_t) * (1 + var_size_expansion_factor) :
                          tiledb::impl::type_size(dim.type());
        } else {
            throw std::runtime_error(fmt::format("[MemoryPoolAllocationStrategy] Missing column name '{}'", column));
        }
    }

    // Ensure minimum buffer size is multiple of 8
    buffer_unit_size = (memory_budget / weight / 8) * 8;
}

std::pair<size_t, size_t> MemoryPoolAllocationStrategy::get_buffer_sizes(
    const std::variant<tiledb::Attribute, tiledb::Dimension> column) const {
    return std::visit(
        [&](auto&& arg) {
            using T = std::decay_t<decltype(arg)>;

            bool is_var = false;
            if constexpr (std::is_same_v<T, tiledb::Attribute>) {
                is_var = arg.variable_sized();
            } else {
                is_var = arg.cell_val_num() == TILEDB_VAR_NUM || arg.type() == TILEDB_STRING_ASCII ||
                         arg.type() == TILEDB_STRING_UTF8;
            }

            size_t column_budget = (is_var ? sizeof(uint64_t) * var_size_expansion_factor :
                                             tiledb::impl::type_size(arg.type())) *
                                   buffer_unit_size;
            size_t num_cells = buffer_unit_size;

            return std::make_pair(column_budget, num_cells);
        },
        column);
}
}  // namespace tiledbsoma