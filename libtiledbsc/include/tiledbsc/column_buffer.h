#ifndef COLUMN_BUFFER_H
#define COLUMN_BUFFER_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <span>
#include <tiledb/tiledb>

#include "tiledbsc/logger_public.h"

namespace tiledbsc {

using namespace tiledb;

/**
 * @brief Class to store data for a TileDB dimension or attribute.
 *
 */
class ColumnBuffer {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a ColumnBuffer from an array and column name.
     *
     * @param array TileDB array
     * @param name TileDB dimension or attribute name
     * @param num_cells Number of cells to allocate
     * @param bytes_per_cell Bytes per cell for variable length data (optional)
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> create(
        std::shared_ptr<Array> array,
        std::string_view name,
        size_t num_cells,
        std::optional<size_t> bytes_per_cell = std::nullopt);

    /**
     * @brief Create a ColumnBuffer from a dimension.
     *
     * @param dim TileDB dimension
     * @param num_cells Number of cells to allocate
     * @param bytes_per_cell Bytes per cell for variable length data (optional)
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> create(
        const Dimension& dim,
        size_t num_cells,
        std::optional<size_t> bytes_per_cell = std::nullopt);

    /**
     * @brief Create a ColumnBuffer from an attribute.
     *
     * @param attr TileDB attribute
     * @param num_cells Number of cells to allocate
     * @param bytes_per_cell Bytes per cell for variable length data (optional)
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> create(
        const Attribute& attr,
        size_t num_cells,
        std::optional<size_t> bytes_per_cell = std::nullopt);

    /**
     * @brief Create a ColumnBuffer from data.
     *
     * @param name Column name
     * @param type TileDB datatype
     * @param num_cells Number of cells
     * @param data View of data
     * @param offsets View of offsets (optional)
     * @param validity View of validity (optional)
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> create(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        std::span<std::byte> data,
        std::span<uint64_t> offsets = {},
        std::span<uint8_t> validity = {});

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new ColumnBuffer object
     *
     * @param name Column name
     * @param type TileDB datatype
     * @param num_cells Number of cells
     * @param data View of data
     * @param offsets View of offsets (optional)
     * @param validity View of validity (optional)
     */
    ColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        std::vector<std::byte> data,
        std::vector<uint64_t> offsets,
        std::vector<uint8_t> validity);

    /**
     * @brief Destroy the ColumnBuffer object
     *
     */
    ~ColumnBuffer();

    /**
     * @brief Attach this ColumnBuffer to a TileDB query.
     *
     * @param query TileDB query
     */
    void attach(Query& query);

    /**
     * @brief Size num_cells_ to match the read query results.
     *
     * @param query TileDB query
     */
    size_t update_size(const Query& query);

    /**
     * @brief Return a view of the ColumnBuffer data.
     *
     * @tparam T Data type
     * @return std::span<T> data view
     */
    template <typename T>
    std::span<T> data() {
        return std::span<T>((T*)data_.data(), num_cells_);
    }

    /**
     * @brief Return data in a vector of strings.
     *
     * @return std::vector<std::string>
     */
    std::vector<std::string> strings();

    /**
     * @brief Return a string_view of the string at the provided cell index.
     *
     * @param index Cell index
     * @return std::string_view string view
     */
    std::string_view string_view(uint64_t index);

    /**
     * @brief Return a view of the ColumnBuffer offsets.
     *
     * @return std::span<uint64_t> offsets view
     */
    std::span<uint64_t> offsets() {
        return std::span<uint64_t>(offsets_);
    }

    /**
     * @brief Return a view of the validity buffer.
     *
     * @return std::span<uint8_t> validity view
     */
    std::span<uint8_t> validity() {
        return std::span<uint8_t>(validity_);
    }

    /**
     * @brief Return the name of the buffer.
     *
     * @return std::string name
     */
    std::string name() {
        return name_;
    }

    /**
     * @brief Return true if the buffer contains variable length data.
     */
    bool is_var() {
        return !offsets_.empty();
    }

    /**
     * @brief Return true if the buffer contains nullable data.
     */
    bool is_nullable() {
        return !validity_.empty();
    }

   private:
    //===================================================================
    //= private static
    //===================================================================

    /**
     * @brief Allocate and return a ColumnBuffer.
     *
     * @param name Column name
     * @param type TileDB datatype
     * @param num_cells Number of cells
     * @param is_var True if variable length data
     * @param is_nullable True if nullable data
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> alloc(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        bool is_var,
        bool is_nullable,
        std::optional<size_t> bytes_per_cell = std::nullopt);

    //===================================================================
    //= private non-static
    //===================================================================

    // Name of the column from the schema.
    std::string name_;

    // Data type of the column from the schema.
    tiledb_datatype_t type_;

    // Bytes per element.
    uint64_t type_size_;

    // Number of cells.
    uint64_t num_cells_;

    // Data buffer.
    std::vector<std::byte> data_;

    // Offsets buffer (optional).
    std::vector<uint64_t> offsets_;

    // Validity buffer (optional).
    std::vector<uint8_t> validity_;
};

}  // namespace tiledbsc
#endif
