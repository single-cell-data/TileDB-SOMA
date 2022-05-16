#ifndef TILEDBSC_BUFFER_SET_H
#define TILEDBSC_BUFFER_SET_H

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <tiledb/tiledb>
#include "tiledbsc_export.h"

namespace tiledbsc {

using DELEM_T = std::byte;

/**
 *
 * Struct representing a set of buffers holding the query result data
 * of a TileDB attribute or dimension.
 *
 * @details
 * The object contains a set of buffers to hold data from a single
 * attribute/column. Offsets and validity buffers are optional and
 * may not be initialized depending on the schema.
 */
struct BufferSet {
    BufferSet() = delete;

    /**
     * TODO
     */
    BufferSet(
        std::string name,
        size_t initial_allocation,
        size_t elem_nbytes,
        bool isvar,
        bool isnullable);

    static std::shared_ptr<BufferSet> from_attribute(
        const tiledb::Attribute& attr, size_t nelem);

    static std::shared_ptr<BufferSet> from_dimension(
        const tiledb::Dimension& attr, size_t nelem);

    /**
     * Resize the buffer set hold `new_size` data *elements*
     */
    void resize(size_t new_size);

    /**
     * Returns true if the buffer set represents a variable-length field
     */
    bool isvar();

    /**
     * Returns true if the buffer set represents a nullable field
     */
    bool isnullable();

    /**
     * Returns number of bytes for each element in the buffer
     */
    size_t elem_nbytes();

    /**
     * Returns name of attribute/dimension represented by this bufferset
     */
    std::string name();

    /* internal methods */
    size_t nelem_data_nbytes(size_t nelem);
    size_t nelem_offsets_nbytes(size_t nelem);
    size_t nelem_validity_nbytes(size_t nelem);

    /** properties **/
    std::string name_;
    size_t elem_nbytes_;

    /** data members **/
    std::vector<DELEM_T> data;
    std::optional<std::vector<DELEM_T>> offsets;
    std::optional<std::vector<DELEM_T>> validity;
};

};  // namespace tiledbsc

#endif