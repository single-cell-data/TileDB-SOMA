#include "tiledbsc/column_buffer.h"
#include "tiledbsc/logger_public.h"

namespace tiledbsc {

using namespace tiledb;

ColumnBuffer ColumnBuffer::create(
    const Array& array,
    std::string_view name,
    size_t num_cells,
    std::optional<size_t> bytes_per_cell) {
    auto name_str = std::string(name);  // string for TileDB API
    auto schema = array.schema();

    if (schema.has_attribute(name_str)) {
        return create(schema.attribute(name_str), num_cells, bytes_per_cell);
    } else if (schema.domain().has_dimension(name_str)) {
        return create(
            schema.domain().dimension(name_str), num_cells, bytes_per_cell);
    }

    throw TileDBSCError("[ColumnBuffer] Column name not found: " + name_str);
}

ColumnBuffer ColumnBuffer::create(
    const tiledb::Dimension& dim,
    size_t num_cells,
    std::optional<size_t> bytes_per_cell) {
    bool is_var = dim.cell_val_num() == TILEDB_VAR_NUM ||
                  dim.type() == TILEDB_STRING_ASCII ||
                  dim.type() == TILEDB_STRING_UTF8;

    return ColumnBuffer::alloc(
        dim.name(), dim.type(), num_cells, is_var, false, bytes_per_cell);
}

ColumnBuffer ColumnBuffer::create(
    const tiledb::Attribute& attr,
    size_t num_cells,
    std::optional<size_t> bytes_per_cell) {
    bool is_var = attr.cell_val_num() == TILEDB_VAR_NUM;

    return ColumnBuffer::alloc(
        attr.name(),
        attr.type(),
        num_cells,
        is_var,
        attr.nullable(),
        bytes_per_cell);
}

ColumnBuffer ColumnBuffer::create(
    std::string_view name,
    tiledb_datatype_t type,
    std::span<std::byte> data,
    std::span<uint64_t> offsets,
    std::span<uint8_t> validity) {
    if (!offsets.empty() && data.size() != offsets.size()) {
        // TODO
    }
    if (!validity.empty() && data.size() != validity.size()) {
        // TODO
    }

    std::vector<std::byte> d{data.begin(), data.end()};
    std::vector<uint64_t> o{offsets.begin(), offsets.end()};
    std::vector<uint8_t> v{validity.begin(), validity.end()};

    return ColumnBuffer(name, type, d, o, v);
}

ColumnBuffer ColumnBuffer::alloc(
    std::string_view name,
    tiledb_datatype_t type,
    size_t num_cells,
    bool is_var,
    bool is_nullable,
    std::optional<size_t> bytes_per_cell) {
    // Compute the number of bytes for the data buffer using `bytes_per_cell` if
    // provided, otherwise use the type_size.
    auto type_size = tiledb::impl::type_size(type);
    auto num_bytes = num_cells * (bytes_per_cell ? *bytes_per_cell : type_size);
    auto data = std::vector<std::byte>(num_bytes);

    // TileDB requires num_cells offset values. We allocate num_cells + 1 to
    // match the Arrow offsets format and enable conversion between the TileDB
    // and Arrow formats.
    // https://arrow.apache.org/docs/format/Columnar.html#variable-size-binary-layout
    auto offsets = is_var ? std::vector<uint64_t>(num_cells + 1) :
                            std::vector<uint64_t>{};

    auto validity = is_nullable ? std::vector<uint8_t>(num_cells) :
                                  std::vector<uint8_t>{};

    return ColumnBuffer(name, type, data, offsets, validity);
}

ColumnBuffer::ColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    std::vector<std::byte> data,
    std::vector<uint64_t> offsets,
    std::vector<uint8_t> validity)
    : name_(name)
    , type_(type)
    , type_size_(tiledb::impl::type_size(type))
    , data_(data)
    , offsets_(offsets)
    , validity_(validity) {
}

ColumnBuffer::~ColumnBuffer() = default;

}  // namespace tiledbsc
