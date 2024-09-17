#ifndef ARROW_ADAPTER_H
#define ARROW_ADAPTER_H

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

// https://arrow.apache.org/docs/format/CDataInterface.html
// https://arrow.apache.org/docs/format/Columnar.html#buffer-listing-for-each-layout
// https://arrow.apache.org/docs/format/CDataInterface.html#exporting-a-simple-int32-array

#include "nanoarrow/nanoarrow.hpp"
#include "nlohmann/json.hpp"

namespace tiledbsoma {

using namespace tiledb;
using json = nlohmann::json;

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
    /* Set the ZstdFilter's level for DataFrame dims */
    int32_t dataframe_dim_zstd_level = 3;

    /* Set the ZstdFilter's level for SparseNDArray dims */
    int32_t sparse_nd_array_dim_zstd_level = 3;

    /* Set the ZstdFilter's level for DenseNDArray dims */
    int32_t dense_nd_array_dim_zstd_level = 3;

    /* Set whether to write the X data chunked */
    bool write_X_chunked = true;

    /* Set the goal chunk nnz */
    uint64_t goal_chunk_nnz = 100000000;

    /* Server-side parameter to set the cap nbytes */
    uint64_t remote_cap_nbytes = 2400000000;

    /* Set the tile capcity for sparse arrays */
    uint64_t capacity = 100000;

    /**
     * Available filters with associated options are
     * [
     *     {
     *         "name": "GZIP", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "ZSTD", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "LZ4", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "BZIP2", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "RLE", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "DELTA",
     *         "COMPRESSION_LEVEL": (int32_t),
     *         "COMPRESSION_REINTERPRET_DATATYPE": (uint8_t)
     *     },
     *     {
     *         "name": "DOUBLE_DELTA",
     *         "COMPRESSION_LEVEL": (int32_t),
     *         "COMPRESSION_REINTERPRET_DATATYPE": (uint8_t)
     *     },
     *     {
     *         "name": "BIT_WIDTH_REDUCTION",
     *         "BIT_WIDTH_MAX_WINDOW": (uint32_t)
     *     },
     *     {
     *         "name": "POSITIVE_DELTA", "POSITIVE_DELTA_MAX_WINDOW":
     * (uint32_t),
     *     },
     *     {
     *         "name": "DICTIONARY_ENCODING", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "SCALE_FLOAT",
     *         "SCALE_FLOAT_FACTOR": (double),
     *         "SCALE_FLOAT_OFFSET": (double),
     *         "SCALE_FLOAT_BYTEWIDTH": (uint64_t),
     *     },
     *     {
     *         "name": "WEBP",
     *         "WEBP_INPUT_FORMAT": (uint8_t),
     *         "WEBP_QUALITY": (float),
     *         "WEBP_LOSSLESS": (uint8_t),
     *     },
     *     "CHECKSUM_MD5",
     *     "CHECKSUM_SHA256",
     *     "XOR",
     *     "BITSHUFFLE",
     *     "BYTESHUFFLE",
     *     "NOOP"
     * ]
     *
     */
    std::string
        offsets_filters = R"(["DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD"])";

    /* Set the validity filters. */
    std::string validity_filters = "";

    /* Set whether the TileDB Array allows duplicate values */
    bool allows_duplicates = false;

    /* Set the tile order as "row", "row-major", "col", or "col-major" */
    std::optional<std::string> tile_order = std::nullopt;

    /* Set the cell order as "hilbert", "row", "row-major", "col", or
     * "col-major"
     */
    std::optional<std::string> cell_order = std::nullopt;

    /* Set the filters for attributes.
     *
     * Example:
     * {
     *     "attr_name": {
     *          "filters": ["XOR", {"name": "GZIP", "COMPRESSION_LEVEL": 3}]
     *     }
     * }
     *
     */
    std::string attrs = "";

    /* Set the filters and tiles for dimensions.
     *
     * Example:
     * {
     *     "dim_name": {"filters": ["NoOpFilter"], "tile": 8}
     * }
     *
     */

    std::string dims = "";

    /* Set whether the array should be consolidated and vacuumed after writing
     */
    bool consolidate_and_vacuum = false;
};

/**
 * This is our application-specific wrapper around nanoarrow.
 *
 * Nominal use-cases include full dataframe-wide schemas from Python/R for
 * DataFrame::create, domain/extent/etc information from Python/R for
 * setting core dimensions, and readback to Python/R of core domain, current
 * domain, and non-empty domain.
 *
 * nanoarrow is crucial to our code-centralization efforts in libtiledbsoma.
 * Python and R can trivially do heterogenous lists like ["abc", 123, true]
 * or list("abc", 123, TRUE); C++ cannot -- at least not so compactly. The
 * combination of ArrowSchema (for the types) and ArrowArray (for the data) is
 * key for Python/R interop with C++ libtiledbsoma.
 */

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
     * The number of rows in index_column_info was three without core
     * current-domain support, and is five with core current-domain support:
     *
     * - Slot 0: core domain low value (inclusive)
     * - Slot 1: core domain high value (inclusive)
     * - Slot 2: core extent parameter
     * - Slot 3: core current-domain low value (inclusive)
     * - Slot 4: core current-domain high value (inclusive)
     *
     * @return tiledb::ArraySchema
     */
    static ArraySchema tiledb_schema_from_arrow_schema(
        std::shared_ptr<Context> ctx,
        std::unique_ptr<ArrowSchema> arrow_schema,
        ArrowTable index_column_info,
        std::string soma_type,
        bool is_sparse = true,
        PlatformConfig platform_config = PlatformConfig());

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

    /**
     * @brief This is a keystroke-saver.
     */
    static std::string tdb_to_arrow_type(tiledb_datatype_t tiledb_datatype) {
        return std::string(to_arrow_format(tiledb_datatype));
    }

    /**
     * @brief Creates a nanoarrow ArrowSchema given the names and TileDB
     * datatypes.
     *
     * Note that the parents and children in nanoarrow are both of type
     * ArrowSchema. This constructs the parent and the children.
     */
    static std::unique_ptr<ArrowSchema> make_arrow_schema(
        const std::vector<std::string>& names,
        const std::vector<tiledb_datatype_t>& tiledb_datatypes);

    /**
     * @brief Creates a nanoarrow ArrowArray which accommodates
     * a varying number of columns.
     *
     * Note that the parents and children in nanoarrow are both of type
     * ArrowArray. This constructs the parent and not the children.
     */
    static std::unique_ptr<ArrowArray> make_arrow_array_parent(int num_columns);

    /**
     * @brief Creates a nanoarrow ArrowArray for a single column.
     *
     * Note that the parents and children in nanoarrow are both of type
     * ArrowArray. This constructs not the parent but rather a child.
     */
    template <typename T>
    static ArrowArray* make_arrow_array_child(const std::pair<T, T>& pair) {
        std::vector<T> v({pair.first, pair.second});
        return make_arrow_array_child<T>(v);
    }

    template <typename T>
    static ArrowArray* make_arrow_array_child(const std::vector<T>& v) {
        // Use new here, not malloc, to match ArrowAdapter::release_array
        auto arrow_array = new ArrowArray;

        size_t n = v.size();

        arrow_array->length = n;
        arrow_array->null_count = 0;
        arrow_array->offset = 0;
        arrow_array->n_buffers = 2;
        arrow_array->release = &ArrowAdapter::release_array;
        arrow_array->buffers = new const void*[2];
        arrow_array->n_children = 0;  // leaf/child node

        arrow_array->buffers[0] = nullptr;
        // Use malloc here, not new, to match ArrowAdapter::release_array
        T* dest = (T*)malloc(n * sizeof(T));
        for (size_t i = 0; i < n; i++) {
            dest[i] = v[i];
        }
        arrow_array->buffers[1] = (void*)dest;

        return arrow_array;
    }

    static ArrowArray* make_arrow_array_child(
        const std::vector<std::string>& v) {
        // Use new here, not malloc, to match ArrowAdapter::release_array
        auto arrow_array = new ArrowArray;

        size_t n = v.size();

        arrow_array->length = n;
        arrow_array->null_count = 0;
        arrow_array->offset = 0;
        arrow_array->n_buffers = 2;
        arrow_array->release = &ArrowAdapter::release_array;
        arrow_array->buffers = new const void*[2];
        arrow_array->n_children = 0;  // leaf/child node

        // For core domain, these are always nullptr for strings and cannot be
        // anything else. More general use of this class is WIP on
        // https://github.com/single-cell-data/TileDB-SOMA/issues/2407
        arrow_array->buffers[0] = nullptr;
        arrow_array->buffers[1] = nullptr;

        return arrow_array;
    }

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

    static Dimension _create_dim(
        tiledb_datatype_t type,
        std::string name,
        const void* buff,
        std::shared_ptr<Context> ctx);

    static void _set_current_domain_slot(
        tiledb_datatype_t type,
        const void* buff,
        NDRectangle& ndrect,
        std::string name);

    template <typename T>
    static Dimension _create_dim_aux(
        std::shared_ptr<Context> ctx, std::string name, T* b) {
        return Dimension::create<T>(*ctx, name, {b[0], b[1]}, b[2]);
    }

    static bool _isvar(const char* format);

    static FilterList _create_filter_list(
        std::string filters, std::shared_ptr<Context> ctx);

    static FilterList _create_filter_list(
        json filters, std::shared_ptr<Context> ctx);

    static FilterList _create_attr_filter_list(
        std::string name,
        PlatformConfig platform_config,
        std::shared_ptr<Context> ctx);

    static FilterList _create_dim_filter_list(
        std::string name,
        PlatformConfig platform_config,
        std::string soma_type,
        std::shared_ptr<Context> ctx);

    static Filter _get_zstd_default(
        PlatformConfig platform_config,
        std::string soma_type,
        std::shared_ptr<Context> ctx);

    static void _append_to_filter_list(
        FilterList filter_list, json filter, std::shared_ptr<Context> ctx);

    static void _set_filter_option(
        Filter filter, std::string option_name, json value);

    static tiledb_layout_t _get_order(std::string order);

};  // class ArrowAdapter

};  // namespace tiledbsoma
#endif
