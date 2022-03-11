#include <memory>
#include <vector>
#include <optional>

// TODO remove work-around in current libtiledb
#include <stdexcept>
#include <tiledb/tiledb>

#include <tiledbsc/buffer_set.h>

namespace {
    bool isvar(const tiledb::Attribute& attr) {
        return attr.cell_val_num() == TILEDB_VAR_NUM;
    }
    bool isvar(const tiledb::Dimension& dim) {
        return (dim.cell_val_num() == TILEDB_VAR_NUM ||
                dim.type() == TILEDB_STRING_ASCII ||
                dim.type() == TILEDB_STRING_UTF8);
    }
}

namespace tiledbsc {

/* ********************************* */
/*              BufferSet            */
/* ********************************* */

BufferSet::BufferSet(
    std::string name,
    size_t nelem, size_t elem_nbytes,
    bool isvar, bool isnullable)
    :
    name_(name), elem_nbytes_(elem_nbytes), offsets(std::nullopt), validity(std::nullopt)
{
    data = std::vector<DELEM_T>(nelem * elem_nbytes);

    if (isvar) {
        offsets = std::vector<DELEM_T>(nelem_offsets_nbytes(nelem));
    }
    if (isnullable) {
        validity = std::vector<DELEM_T>(nelem_validity_nbytes(nelem));
    }
}

size_t BufferSet::nelem_data_nbytes(size_t nelem) {
    return nelem * elem_nbytes_;
}

size_t BufferSet::nelem_offsets_nbytes(size_t nelem) {
    return nelem * sizeof(uint64_t);
}

size_t BufferSet::nelem_validity_nbytes(size_t nelem) {
    return nelem * sizeof(uint8_t);
}

bool BufferSet::isvar() {
    return offsets.has_value();
}

bool BufferSet::isnullable() {
    return validity.has_value();
}

size_t BufferSet::elem_nbytes() {
    return elem_nbytes_;
}

std::string BufferSet::name() {
    return name_;
}

void BufferSet::resize(size_t nelem) {
    data.resize(nelem_data_nbytes(nelem));

    if (offsets) {
        offsets.value().resize(nelem_offsets_nbytes(nelem));
    }

    if (validity) {
        validity.value().resize(nelem_offsets_nbytes(nelem));
    }
}

std::shared_ptr<BufferSet> BufferSet::from_attribute(
    const tiledb::Attribute& attr, size_t nelem
)
{
    size_t elem_nbytes = tiledb::impl::type_size(attr.type());

    return std::make_shared<BufferSet>(
        attr.name(),
        nelem,
        elem_nbytes,
        ::isvar(attr),
        attr.nullable()
    );
}

std::shared_ptr<BufferSet> BufferSet::from_dimension(
    const tiledb::Dimension& dim, size_t nelem
)
{
    size_t elem_nbytes = tiledb::impl::type_size(dim.type());

    return std::make_unique<BufferSet>(
        dim.name(),
        nelem,
        elem_nbytes,
        ::isvar(dim),
        false
    );
}

}; // namespace tiledb