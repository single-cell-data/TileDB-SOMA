#ifndef TILEDBSC_BUFFER_SET_H
#define TILEDBSC_BUFFER_SET_H

#include "tiledbsc_export.h"

#include <vector>
#include <string>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>

#include <tiledb/tiledb>

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
class TILEDBSC_EXPORT BufferSet {
public:
    BufferSet(
        const std::string& name,
        tiledb_datatype_t datatype,
        size_t elem_nbytes,
        std::vector<DELEM_T> data,
        std::optional<std::vector<uint64_t>> offsets,
        std::optional<std::vector<DELEM_T>> validity,
        bool convert_bitmap = false
    );

    ~BufferSet();

    static BufferSet from_attribute(
        const tiledb::Attribute& attr, size_t data_nelem
    );

    static BufferSet from_dimension(
        const tiledb::Dimension& attr, size_t data_nelem
    );

    /**
     * Construct a BufferSet from data.
     * Create as nullable if the optional
     * validity buffer is supplied.
     */
    static BufferSet from_data(
        const std::string& name,
        tiledb_datatype_t type,
        size_t elem_nbytes,
        std::span<std::byte> data,
        std::optional<std::span<uint64_t>> offsets = std::nullopt,
        std::optional<std::span<std::byte>> validity = std::nullopt
    );

    /**
     * Construct desired BufferSet with initial element count.
     */
    static BufferSet alloc(
        const std::string& name,
        tiledb_datatype_t type,
        size_t initial_nelem,
        size_t elem_nbytes,
        bool isvar,
        bool isnullable);

    /**
     * Return view of data
     */
    template <typename T>
    std::span<T> data() {
        return std::span<T>(data_);
    }

    /**
     * Return view of data
     */
    std::span<uint64_t> offsets();

    /**
     * Return view of validity bytemap
     */
    std::span<std::byte> validity();

    /**
     * Resize the buffer set hold `new_size` data *elements*
     * and optionally `new_ncells` *cells* (offsets, validity)
     */
    void resize(size_t new_data_nelem, std::optional<size_t> new_ncells = std::nullopt);

    /**
     * Returns true if the buffer set represents a variable-length field
     */
    // TODO rename this to isvarlength?
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
     * Returns number of cells in the bufferset
     */
    size_t num_cells();

    /**
     * Returns number of null values in the bufferset
     */
    size_t num_nulls();

    /**
     * Sets the number of null values in the bufferset
     */
    void set_num_nulls(size_t);

    /**
     * Returns name of attribute/dimension represented by this bufferset
     */
    std::string name();

    /**
     * Returns the datatype of the buffer
     */
    tiledb_datatype_t datatype();

    /* internal methods */
    static size_t nelem_data_nbytes(size_t nelem);
    static size_t nelem_offsets_nbytes(size_t nelem);
    static size_t nelem_validity_nbytes(size_t nelem);

    /** properties **/
    std::string name_;
    tiledb_datatype_t datatype_;
    size_t elem_nbytes_;
    size_t num_nulls_;

    /** data members **/
    std::vector<DELEM_T> data_;
    std::optional<std::vector<uint64_t>> offsets_;
    std::optional<std::vector<DELEM_T>> validity_;
};

}; // namespace tiledb

#endif