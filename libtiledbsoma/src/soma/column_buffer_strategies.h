#ifndef COLUMN_BUFFER_STRATEGIES_H
#define COLUMN_BUFFER_STRATEGIES_H

#include <tiledb/tiledb>

#include <span>
#include <utility>
#include <variant>

namespace tiledbsoma {

inline constexpr std::string_view CONFIG_KEY_INIT_BYTES = "soma.init_buffer_bytes";

inline constexpr std::string_view CONFIG_KEY_MEMORY_POOL_SIZE{"soma.read.memory_budget"};
inline constexpr std::string_view CONFIG_KEY_VAR_SIZED_FACTOR{"soma.read.var_size_factor"};

class ColumnBufferAllocationStrategy {
   public:
    virtual ~ColumnBufferAllocationStrategy() = default;

    virtual std::pair<size_t, size_t> get_buffer_sizes(
        std::variant<tiledb::Attribute, tiledb::Dimension> column) const = 0;
};

class BasicAllocationStrategy : public ColumnBufferAllocationStrategy {
   public:
    BasicAllocationStrategy(const tiledb::Array& array);

    std::pair<size_t, size_t> get_buffer_sizes(
        std::variant<tiledb::Attribute, tiledb::Dimension> column) const override;

   private:
    size_t buffer_size = 1 << 30;
};

inline constexpr size_t DEFAULT_MEMORY_POOL = 1 << 30;
inline constexpr size_t DEFAULT_VAR_SIZE_FACTOR = 2;

class MemoryPoolAllocationStrategy : public ColumnBufferAllocationStrategy {
   public:
    MemoryPoolAllocationStrategy(std::span<std::string> columns, const tiledb::Array& array);

    std::pair<size_t, size_t> get_buffer_sizes(
        std::variant<tiledb::Attribute, tiledb::Dimension> column) const override;

   private:
    size_t buffer_unit_size = 1 << 30;
    size_t var_size_expansion_factor = 2;
};
}  // namespace tiledbsoma
#endif