#include <memory>
#include <optional>
#include <vector>

// TODO remove work-around in current libtiledb
#include <stdexcept>
#include <tiledb/tiledb>

#include <tiledbsc/buffer_set.h>
#include <tiledbsc/common.h>
#include <tiledbsc/util.h>

namespace {
bool isvar(const tiledb::Attribute& attr) {
    return attr.cell_val_num() == TILEDB_VAR_NUM;
}
bool isvar(const tiledb::Dimension& dim) {
    return (
        dim.cell_val_num() == TILEDB_VAR_NUM ||
        dim.type() == TILEDB_STRING_ASCII || dim.type() == TILEDB_STRING_UTF8);
}
}  // namespace

namespace tiledbsc {

using namespace std;

/* ********************************* */
/*              BufferSet            */
/* ********************************* */

BufferSet::BufferSet(
    const std::string& name,
    tiledb_datatype_t datatype,
    size_t elem_nbytes,
    std::vector<DELEM_T> data,
    std::optional<std::vector<uint64_t>> offsets,
    std::optional<std::vector<DELEM_T>> validity,
    bool convert_bitmap)
    : name_(name)
    , datatype_(datatype)
    , elem_nbytes_(elem_nbytes)
    , data_(data)
    , offsets_(offsets)
    , validity_(validity) {
    if (elem_nbytes == 0) {
        throw TileDBSCError("[BufferSet] Invalid elem_nbytes in constructor");
    }

    if (tiledb_datatype_size(datatype) != elem_nbytes) {
        throw TileDBSCError(
            "[BufferSet] elem_nbytes does not match datatype size");
    }
    if (offsets && validity) {
        if (offsets.value().size() - 1 != validity.value().size())
            throw TileDBSCError(
                "[BufferSet] Mismatched validity and offsets buffer size");
    };
    if (!offsets && validity) {
        size_t data_ncells = data.size() / elem_nbytes;
        if (data_ncells != validity.value().size())
            throw TileDBSCError(
                "[BufferSet] Mismatched data and validity buffer size");
    }

    if (convert_bitmap) {
        // util::
    }
}

BufferSet::~BufferSet() = default;

size_t BufferSet::nelem_offsets_nbytes(size_t nelem) {
    return nelem * sizeof(uint64_t);
}

size_t BufferSet::nelem_validity_nbytes(size_t nelem) {
    return nelem * sizeof(uint8_t);
}

bool BufferSet::isvar() {
    return offsets_.has_value();
}

bool BufferSet::isnullable() {
    return validity_.has_value();
}

size_t BufferSet::elem_nbytes() {
    return elem_nbytes_;
}

size_t BufferSet::num_cells() {
    if (isvar())
        // always Arrow offsets
        return offsets_.value().size() - 1;
    else
        return data_.size() / elem_nbytes_;
}

void BufferSet::set_num_nulls(size_t nnull) {
    if (!isnullable())
        throw std::runtime_error(
            "num_nulls is only valid for nullable buffersets");

    num_nulls_ = nnull;
}

size_t BufferSet::num_nulls() {
    if (!isnullable())
        throw std::runtime_error(
            "num_nulls is only valid for nullable buffersets");

    return num_nulls_;
}

tiledb_datatype_t BufferSet::datatype() {
    return datatype_;
}

std::string BufferSet::name() {
    return name_;
}

std::span<uint64_t> BufferSet::offsets() {
    if (!offsets_.has_value())
        throw std::logic_error(
            "Offsets buffer not defined for this BufferSet!");

    return span<uint64_t>(offsets_.value());
}

span<byte> BufferSet::validity() {
    if (!validity_.has_value())
        throw std::logic_error(
            "Validity buffer not defined for this BufferSet!");

    return validity_.value();
}

void BufferSet::resize(size_t data_nelem, std::optional<size_t> ncells) {
    /* TODO CHECK NEW SIZES */
    data_.resize(data_nelem * elem_nbytes_);

    if (ncells) {
        if (offsets_) {
            offsets_.value().resize(ncells.value() + 1);
        }
        if (validity_) {
            validity_.value().resize(ncells.value());
        }
    }
}

BufferSet BufferSet::from_attribute(
    const tiledb::Attribute& attr, size_t data_nelem) {
    size_t elem_nbytes = tiledb::impl::type_size(attr.type());

    return BufferSet::alloc(
        attr.name(),
        attr.type(),
        data_nelem,
        elem_nbytes,
        ::isvar(attr),
        attr.nullable());
}

BufferSet BufferSet::from_dimension(
    const tiledb::Dimension& dim, size_t data_nelem) {
    size_t elem_nbytes = tiledb::impl::type_size(dim.type());

    return BufferSet::alloc(
        dim.name(), dim.type(), data_nelem, elem_nbytes, ::isvar(dim), false);
}

BufferSet BufferSet::alloc(
    const std::string& name,
    tiledb_datatype_t datatype,
    size_t data_nelem,
    size_t elem_nbytes,
    bool isvar,
    bool isnullable) {
    auto data = std::vector<DELEM_T>(data_nelem * elem_nbytes);

    auto offsets = isvar ? std::make_optional<std::vector<uint64_t>>(
                               data_nelem + 1) :
                           std::nullopt;

    auto validity = isnullable ?
                        std::make_optional<std::vector<DELEM_T>>(data_nelem) :
                        std::nullopt;

    return BufferSet(name, datatype, elem_nbytes, data, offsets, validity);
}

BufferSet BufferSet::from_data(
    const std::string& name,
    tiledb_datatype_t type,
    size_t elem_nbytes,
    span<byte> data,
    optional<span<uint64_t>> offsets,
    optional<span<byte>> validity) {
    if (offsets)
        if (data.size() > 1 && offsets.value().size() <= 2)
            throw TileDBSCError("Offset count is impossible for given data");

    if (offsets && validity)
        if (offsets.value().size() - 1 != validity.value().size())
            throw std::logic_error(
                "Offsets count does not match validity count");

    vector<byte> d{data.begin(), data.end()};
    auto o = offsets ? make_optional<vector<uint64_t>>(
                           offsets.value().begin(), offsets.value().end()) :
                       nullopt;
    auto v = validity ? make_optional<vector<byte>>(
                            validity.value().begin(), validity.value().end()) :
                        nullopt;

    return BufferSet(name, type, elem_nbytes, d, o, v);
}

};  // namespace tiledbsc