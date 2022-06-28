#include "tiledbsc/column_buffer.h"
#include "tiledbsc/common.h"
#include "tiledbsc/logger_public.h"

namespace tiledbsc {

using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::shared_ptr<ColumnBuffer> ColumnBuffer::create(
    std::shared_ptr<Array> array,
    std::string_view name,
    size_t num_cells,
    std::optional<size_t> bytes_per_cell) {
    auto name_str = std::string(name);  // string for TileDB API
    auto schema = array->schema();

    if (schema.has_attribute(name_str)) {
        return create(schema.attribute(name_str), num_cells, bytes_per_cell);
    } else if (schema.domain().has_dimension(name_str)) {
        return create(
            schema.domain().dimension(name_str), num_cells, bytes_per_cell);
    }

    throw TileDBSCError("[ColumnBuffer] Column name not found: " + name_str);
}

std::shared_ptr<ColumnBuffer> ColumnBuffer::create(
    const tiledb::Dimension& dim,
    size_t num_cells,
    std::optional<size_t> bytes_per_cell) {
    bool is_var = dim.cell_val_num() == TILEDB_VAR_NUM ||
                  dim.type() == TILEDB_STRING_ASCII ||
                  dim.type() == TILEDB_STRING_UTF8;

    return ColumnBuffer::alloc(
        dim.name(), dim.type(), num_cells, is_var, false, bytes_per_cell);
}

std::shared_ptr<ColumnBuffer> ColumnBuffer::create(
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

std::shared_ptr<ColumnBuffer> ColumnBuffer::create(
    std::string_view name,
    tiledb_datatype_t type,
    size_t num_cells,
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

    return std::make_shared<ColumnBuffer>(name, type, num_cells, d, o, v);
}

//===================================================================
//= public non-static
//===================================================================

ColumnBuffer::ColumnBuffer(
    std::string_view name,
    tiledb_datatype_t type,
    size_t num_cells,
    std::vector<std::byte> data,
    std::vector<uint64_t> offsets,
    std::vector<uint8_t> validity)
    : name_(name)
    , type_(type)
    , type_size_(tiledb::impl::type_size(type))
    , num_cells_(num_cells)
    , data_(data)
    , offsets_(offsets)
    , validity_(validity) {
}

ColumnBuffer::~ColumnBuffer(){};

void ColumnBuffer::attach(Query& query) {
    LOG_DEBUG(fmt::format("Attaching buffer {} to query", name_));
    query.set_data_buffer(name_, data_);
    if (!offsets_.empty()) {
        query.set_offsets_buffer(name_, offsets_);
    }
    if (!validity_.empty()) {
        query.set_validity_buffer(name_, validity_);
    }
}

size_t ColumnBuffer::update_size(const Query& query) {
    auto [num_offsets, num_elements] = query.result_buffer_elements()[name_];

    if (is_var()) {
        num_cells_ = num_offsets;
        // Add extra offset for arrow. Resize the offsets buffer if needed.
        if (offsets_.size() < num_offsets + 1) {
            offsets_.resize(num_offsets + 1);
        }
        offsets_[num_offsets] = num_elements;
    } else {
        num_cells_ = num_elements / type_size_;
    }

    return num_cells_;
}

std::vector<std::string> ColumnBuffer::strings() {
    std::vector<std::string> result;

    for (size_t i = 0; i < num_cells_; i++) {
        result.push_back(std::string(string_view(i)));
    }

    return result;
}

std::string_view ColumnBuffer::string_view(uint64_t index) {
    auto start = offsets_[index];
    auto len = offsets_[index + 1] - offsets_[index];
    return std::string_view((char*)(data_.data() + start), len);
}

//===================================================================
//= private static
//===================================================================

std::shared_ptr<ColumnBuffer> ColumnBuffer::alloc(
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

    return std::make_shared<ColumnBuffer>(
        name, type, num_cells, data, offsets, validity);
}

}  // namespace tiledbsc
