#ifndef COLUMN_BUFFER_H
#define COLUMN_BUFFER_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <span>
#include <tiledb/tiledb>
#include <tiledbsc/tiledbsc>

namespace tiledbsc {

using namespace tiledb;

/**
 * @brief Buffer to store data from a TileDB dimension or attribute.
 *
 */
class ColumnBuffer {
   public:
    /**
     * @brief Create a ColumnBuffer from an array and column name.
     *
     * @param array TileDB array
     * @param name TileDB dimension or attribute name
     * @param num_cells Number of cells to allocate
     * @param bytes_per_cell Bytes per cell for variable length data (optional)
     * @return ColumnBuffer
     */
    static ColumnBuffer create(
        const Array& array,
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
    static ColumnBuffer create(
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
    static ColumnBuffer create(
        const Attribute& attr,
        size_t num_cells,
        std::optional<size_t> bytes_per_cell = std::nullopt);

    /**
     * @brief Create a ColumnBuffer from data.
     *
     * @param name Column name
     * @param type TileDB datatype
     * @param data View of data
     * @param offsets View of offsets (optional)
     * @param validity View of validity (optional)
     * @return ColumnBuffer
     */
    static ColumnBuffer create(
        std::string_view name,
        tiledb_datatype_t type,
        std::span<std::byte> data,
        std::span<uint64_t> offsets = {},
        std::span<uint8_t> validity = {});

    /**
     * @brief Attach this ColumnBuffer to a TileDB query.
     *
     * @param query TileDB query
     */
    void attach(Query& query) {
        query.set_data_buffer(name_, data_);
        if (!offsets_.empty()) {
            query.set_offsets_buffer(name_, offsets_);
        }
        if (!validity_.empty()) {
            query.set_validity_buffer(name_, validity_);
        }
    }

    /**
     * @brief Resize the buffers to match the read query results.
     *
     * @param query TileDB query
     */
    void resize(const Query& query) {
        auto [num_offsets, num_elements] = query
                                               .result_buffer_elements()[name_];

        // Resize data buffer (bytes) to match query result size (elements)
        data_.resize(num_elements * type_size_);

        if (num_offsets) {
            // Variable length data
            offsets_.resize(num_offsets);
            // Extra offset value for arrow
            offsets_.push_back(data_.size());
            validity_.resize(num_offsets);
        } else {
            // Fixed length data
            validity_.resize(num_elements);
        }
    }

    /**
     * @brief Return a view of the ColumnBuffer data.
     *
     * @tparam T Data type
     * @return std::span<T> data view
     */
    template <typename T>
    std::span<T> data() {
        return std::span<T>((T*)data_.data(), data_.size() / type_size_);
    }

    /**
     * @brief Return a string_view of the string at the provided cell index.
     *
     * @param index Cell index
     * @return std::string_view string view
     */
    std::string_view string_view(uint64_t index) {
        auto start = offsets_[index];
        auto len = offsets_[index + 1] - offsets_[index];
        return std::string_view((char*)(data_.data() + start), len);
    }

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
     *
     * @return true
     * @return false
     */
    bool is_var() {
        return !offsets_.empty();
    }

    /**
     * @brief Return true if the buffer contains nullable data.
     *
     * @return true
     * @return false
     */
    bool is_nullable() {
        return !validity_.empty();
    }

    /**
     * @brief Destroy the ColumnBuffer object
     *
     */
    ~ColumnBuffer();

    ColumnBuffer() = delete;

   private:
    /**
     * @brief
     *
     * @param name Column name
     * @param type TileDB datatype
     * @param num_cells Number of cells
     * @param is_var True if variable length data
     * @param is_nullable True if nullable data
     * @return ColumnBuffer
     */
    static ColumnBuffer alloc(
        std::string_view name,
        tiledb_datatype_t type,
        size_t num_cells,
        bool is_var,
        bool is_nullable,
        std::optional<size_t> bytes_per_cell = std::nullopt);

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
        std::vector<std::byte> data,
        std::vector<uint64_t> offsets,
        std::vector<uint8_t> validity);

    // Name of the column from the schema.
    std::string name_;

    // Data type of the column from the schema.
    tiledb_datatype_t type_;

    uint64_t type_size_;

    // Data buffer.
    std::vector<std::byte> data_;

    // Offsets buffer (optional).
    std::vector<uint64_t> offsets_;

    // Validity buffer (optional).
    std::vector<uint8_t> validity_;
};

}  // namespace tiledbsc
#endif
