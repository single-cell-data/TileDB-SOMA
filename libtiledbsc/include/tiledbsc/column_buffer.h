#ifndef COLUMN_BUFFER_H
#define COLUMN_BUFFER_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <span>
#include <tiledb/tiledb>

#include "tiledbsc/common.h"
#include "tiledbsc/logger_public.h"

namespace tiledbsc {

using namespace tiledb;

class ColumnBuffer;  // forward declaration

// Using map instead of unordered map to preserve lexigraphical ordering.
// However, we may want to preserve the order of insertion instead.

// Map: column name -> ColumnBuffer
using ArrayBuffers = std::map<std::string, std::shared_ptr<ColumnBuffer>>;

// Map: array name -> ArrayBuffers
using MultiArrayBuffers = std::map<std::string, ArrayBuffers>;

/**
 * @brief Class to store data for a TileDB dimension or attribute.
 *
 */
class ColumnBuffer {
    inline static const size_t DEFAULT_ALLOC_BYTES = 1 << 24;  // 16 MiB
    inline static const std::string
        CONFIG_KEY_INIT_BYTES = "soma.init_buffer_bytes";

   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a ColumnBuffer from an array and column name.
     *
     * @param array TileDB array
     * @param name TileDB dimension or attribute name
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> create(
        std::shared_ptr<Array> array, std::string_view name);

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
     * @param is_var Column type is variable length
     * @param is_nullable Column can contain null values
     */
    ColumnBuffer(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        size_t data,
        bool is_var = false,
        bool is_nullable = false);

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
     * @brief Return the number of cells in the buffer.
     *
     * @return size_t
     */
    size_t size() const {
        return num_cells_;
    }

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
        if (!is_var_) {
            throw TileDBSCError(
                "[ColumnBuffer] Offsets buffer not defined for " + name_);
        }

        return std::span<uint64_t>(offsets_.data(), num_cells_);
    }

    /**
     * @brief Return a view of the validity buffer.
     *
     * @return std::span<uint8_t> validity view
     */
    std::span<uint8_t> validity() {
        if (!is_nullable_) {
            throw TileDBSCError(
                "[ColumnBuffer] Validity buffer not defined for " + name_);
        }
        return std::span<uint8_t>(validity_.data(), num_cells_);
    }

    /**
     * @brief Return the name of the buffer.
     *
     * @return std::string name
     */
    std::string name() const {
        return name_;
    }

    /**
     * @brief Return the type of the buffer.
     *
     * @return tiledb_datatype_t type
     */
    tiledb_datatype_t type() const {
        return type_;
    }

    /**
     * @brief Return true if the buffer contains variable length data.
     */
    bool is_var() const {
        return is_var_;
    }

    /**
     * @brief Return true if the buffer contains nullable data.
     */
    bool is_nullable() const {
        return is_nullable_;
    }

   private:
    //===================================================================
    //= private static
    //===================================================================

    /**
     * @brief Allocate and return a ColumnBuffer.
     *
     * @param array TileDB array
     * @param name Column name
     * @param type TileDB datatype
     * @param is_var True if variable length data
     * @param is_nullable True if nullable data
     * @return ColumnBuffer
     */
    static std::shared_ptr<ColumnBuffer> alloc(
        std::shared_ptr<Array> array,
        std::string_view name,
        tiledb_datatype_t type,
        bool is_var,
        bool is_nullable);

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

    // If true, the data type is variable length
    bool is_var_;

    // If true, the data is nullable
    bool is_nullable_;

    // Data buffer.
    std::vector<std::byte> data_;

    // Offsets buffer (optional).
    std::vector<uint64_t> offsets_;

    // Validity buffer (optional).
    std::vector<uint8_t> validity_;
};

}  // namespace tiledbsc
#endif
