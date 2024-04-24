#ifndef ARROW_ADAPTER_H
#define ARROW_ADAPTER_H

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

// https://arrow.apache.org/docs/format/CDataInterface.html
// https://arrow.apache.org/docs/format/Columnar.html#buffer-listing-for-each-layout
// https://arrow.apache.org/docs/format/CDataInterface.html#exporting-a-simple-int32-array

#include "nanoarrow/nanoarrow.hpp"

namespace tiledbsoma {

using namespace tiledb;

class ColumnBuffer;

/**
 * @brief The ArrowBuffer holds a shared pointer to a ColumnBuffer, which
 * manages the lifetime of a ColumnBuffer used to back an Arrow array.
 *
 * The ArrowArray.release callback will delete the ArrowBuffer, and
 * automatically decrement the use count of the ColumnBuffer's shared pointer.
 *
 */
struct ArrowBuffer {
    ArrowBuffer(std::shared_ptr<ColumnBuffer> buffer)
        : buffer_(buffer){};

    std::shared_ptr<ColumnBuffer> buffer_;
};

using ArrowTable =
    std::pair<std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>>;

class PlatformConfig {
   public:
    uint64_t dataframe_dim_zstd_level = 3;
    uint64_t sparse_nd_array_dim_zstd_level = 3;
    bool write_X_chunked = true;
    uint64_t goal_chunk_nnz = 100000000;
    uint64_t remote_cap_nbytes = 2400000000;
    uint64_t capacity = 100000;
    std::vector<std::string> offsets_filters = {
        "DoubleDeltaFilter", "BitWidthReductionFilter", "ZstdFilter"};
    std::vector<std::string> validity_filters;
    bool allows_duplicates = false;
    std::optional<std::string> tile_order = std::nullopt;
    std::optional<std::string> cell_order = std::nullopt;
    bool consolidate_and_vacuum = false;
};

class ArrowAdapter {
   public:
    static void release_schema(struct ArrowSchema* schema);
    static void release_array(struct ArrowArray* array);

    static bool _isstr(const char* format);

    /**
     * @brief Convert ColumnBuffer to an Arrow array.
     *
     * @return std::pair<std::unique_ptr<ArrowArray>,
     * std::unique_ptr<ArrowSchema>>
     */
    static std::pair<std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>>
    to_arrow(std::shared_ptr<ColumnBuffer> column);

    /**
     * @brief Create a an ArrowSchema from TileDB Schema
     *
     * @return ArrowSchema
     */
    static std::unique_ptr<ArrowSchema> arrow_schema_from_tiledb_array(
        std::shared_ptr<Context> ctx, std::shared_ptr<Array> tiledb_array);

    /**
     * @brief Create a TileDB ArraySchema from ArrowSchema
     *
     * @return tiledb::ArraySchema
     */
    static ArraySchema tiledb_schema_from_arrow_schema(
        std::shared_ptr<Context> ctx,
        std::unique_ptr<ArrowSchema> arrow_schema,
        ArrowTable index_column_info,
        std::optional<PlatformConfig> platform_config);

    /**
     * @brief Get Arrow format string from TileDB datatype.
     *
     * @param tiledb_dtype TileDB datatype.
     * @return std::string_view Arrow format string.
     */
    static std::string_view to_arrow_format(
        tiledb_datatype_t tiledb_dtype, bool use_large = true);

    /**
     * @brief Get TileDB datatype from Arrow format string.
     *
     * @param datatype TileDB datatype.
     * @return std::string_view Arrow format string.
     */
    static tiledb_datatype_t to_tiledb_format(std::string_view arrow_dtype);

    static enum ArrowType to_nanoarrow_type(std::string_view sv);

   private:
    static std::pair<const void*, std::size_t> _get_data_and_length(
        Enumeration& enmr, const void* dst);

    template <typename T>
    static const void* _fill_data_buffer(std::vector<T> src, const void* dst) {
        auto sz = src.size() * sizeof(T);
        dst = (const void*)malloc(sz);
        std::memcpy((void*)dst, src.data(), sz);
        return dst;
    }

    template <typename T>
    static Dimension _create_dim(Context ctx, std::string name, T* b) {
        return Dimension::create<T>(ctx, name, {b[0], b[1]}, b[2]);
    }

    static bool _isvar(const char* format);
};
};  // namespace tiledbsoma

#endif
