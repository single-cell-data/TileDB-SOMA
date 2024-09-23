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
        ArrowArray* child = make_arrow_array_child<T>(v);
        return child;
    }
    /**
     * Helper for make_arrow_array_child. We need the templating
     * in the header file so it's instantiated at all callsites.
     * But fmt::format logging is hard in header files since
     * #include <logger.h> is relative to the includer, which varies.
     */
    static void log_make_arrow_array_child(ArrowArray* child);

    static ArrowArray* make_arrow_array_child_string(
        const std::pair<std::string, std::string>& pair) {
        std::vector<std::string> v({pair.first, pair.second});
        return make_arrow_array_child_string(v);
    }

    template <typename T>
    static ArrowArray* make_arrow_array_child(const std::vector<T>& v) {
        if (std::is_same_v<T, std::string>) {
            throw std::runtime_error(
                "ArrowAdapter::make_arrow_array_child: template-specialization "
                "failure.");
        }

        // Use malloc here, not new, to match ArrowAdapter::release_array
        auto arrow_array = (ArrowArray*)malloc(sizeof(ArrowArray));

        size_t n = v.size();

        arrow_array->length = n;
        arrow_array->null_count = 0;
        arrow_array->offset = 0;

        // Two-buffer model for non-string data:
        // * Slot 0 is the Arrow validity buffer which we leave null
        // * Slot 1 is data, void* but will be derefrenced as T*
        // * There is no offset information
        arrow_array->n_buffers = 2;
        arrow_array->n_children = 0;  // leaf/child node
        arrow_array->buffers = (const void**)malloc(2 * sizeof(void*));
        arrow_array->children = nullptr;  // leaf/child node
        arrow_array->dictionary = nullptr;

        // The nominal use of these methods as of this writing is for
        // low-volume data such as schema information -- less than a
        // kilobyte total. It's simplest and safest to do data copies,
        // for-loop-wise. If we were to extend usage of these methods
        // to bulk data in the megabyte/gigabyte range, we'd want to
        // look at zero-copy for buffers, with variable approaches
        // to memory management.
        arrow_array->release = &ArrowAdapter::release_array;
        arrow_array->private_data = nullptr;

        arrow_array->buffers[0] = nullptr;
        // Use malloc here, not new, to match ArrowAdapter::release_array
        T* dest = (T*)malloc(n * sizeof(T));
        for (size_t i = 0; i < n; i++) {
            dest[i] = v[i];
        }
        arrow_array->buffers[1] = (void*)dest;

        log_make_arrow_array_child(arrow_array);

        return arrow_array;
    }

    // A nominal use of this method is for reporting core domain, current
    // domain, and non-empty domain back to Python/R.  Meanwhile core string
    // dims must always have domain of (nullptr, nullptr); and they have current
    // domain which must _not_ be nullptr pairs.
    //
    // For the former do we give back a column of length 2 with nulls in it,
    // using Arrow's validity buffers?  Or do we use ("", "") as TileDB-Py does?
    //
    // We choose the latter.
    static ArrowArray* make_arrow_array_child_string(
        const std::vector<std::string>& v) {
        // Use malloc here, not new, to match ArrowAdapter::release_array
        auto arrow_array = (ArrowArray*)malloc(sizeof(ArrowArray));

        size_t n = v.size();

        arrow_array->length = n;  // Number of strings, not number of bytes
        arrow_array->null_count = 0;
        arrow_array->offset = 0;

        // Three-buffer model for string data:
        // * Slot 0 is the Arrow uint8_t* validity buffer
        // * Slot 1 is the Arrow offsets buffer: uint32_t* for Arrow string
        //   or uint64_t* for Arrow large_string
        // * Slot 2 is data, void* but will be derefrenced as T*
        arrow_array->n_buffers = 3;
        arrow_array->buffers = (const void**)malloc(3 * sizeof(void*));
        arrow_array->n_children = 0;      // leaf/child node
        arrow_array->children = nullptr;  // leaf/child node
        arrow_array->dictionary = nullptr;

        arrow_array->release = &ArrowAdapter::release_array;
        arrow_array->private_data = nullptr;

        size_t nbytes = 0;
        for (auto e : v) {
            nbytes += e.length();
        }

        // This function produces arrow large_string, which has 64-bit offsets.
        uint64_t* offsets = (uint64_t*)malloc((n + 1) * sizeof(uint64_t));

        // Data
        char* data = (char*)malloc(nbytes * sizeof(char));
        uint64_t dest_start = 0;

        offsets[0] = dest_start;
        for (size_t i = 0; i < n; i++) {
            const std::string& elem = v[i];
            size_t elem_len = elem.size();

            memcpy(&data[dest_start], elem.c_str(), elem_len);
            dest_start += elem_len;
            offsets[i + 1] = dest_start;
        }

        arrow_array->buffers[0] = nullptr;  // validity
        arrow_array->buffers[1] = offsets;
        arrow_array->buffers[2] = data;

        return arrow_array;
    }

    // These table-column getters are, as of this writing, intended primarily
    // for keystroke-reduction in unit-test cases.

    template <typename T>
    static std::vector<T> get_table_column_by_name(
        const ArrowTable& arrow_table, std::string column_name) {
        int64_t index = _get_column_index_from_name(arrow_table, column_name);
        return get_table_column_by_index<T>(arrow_table, index);
    }

    static std::vector<std::string> get_table_string_column_by_name(
        const ArrowTable& arrow_table, std::string column_name) {
        int64_t index = _get_column_index_from_name(arrow_table, column_name);
        return get_table_string_column_by_index(arrow_table, index);
    }

    template <typename T>
    static std::vector<T> get_table_column_by_index(
        const ArrowTable& arrow_table, int64_t column_index) {
        ArrowArray* arrow_array = arrow_table.first.get();
        ArrowSchema* arrow_schema = arrow_table.second.get();
        _check_shapes(arrow_array, arrow_schema);

        if (std::is_same_v<T, std::string>) {
            throw std::runtime_error(
                "SOMAArray::_core_domain_slot: template-specialization "
                "failure.");
        }

        ArrowArray* child = _get_and_check_column(arrow_table, column_index, 2);

        // For our purposes -- reporting domains, etc. -- we don't use the Arrow
        // validity buffers. If this class needs to be extended someday to
        // support arrow-nulls, we can work on that.
        if (child->buffers[0] != nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_column_by_index: validity buffer "
                "unsupported here");
        }

        const void* vdata = child->buffers[1];
        if (vdata == nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_column_by_index: null data buffer");
        }

        const T* data = (T*)vdata;

        std::vector<T> retval(child->length);
        for (auto i = 0; i < child->length; i++) {
            retval[i] = data[i];
        }
        return retval;
    }

    static std::vector<std::string> get_table_string_column_by_index(
        const ArrowTable& arrow_table, int64_t column_index) {
        ArrowArray* arrow_array = arrow_table.first.get();
        ArrowSchema* arrow_schema = arrow_table.second.get();
        _check_shapes(arrow_array, arrow_schema);

        ArrowArray* child = _get_and_check_column(arrow_table, column_index, 3);

        // For our purposes -- reporting domains, etc. -- we don't use the Arrow
        // validity buffers. If this class needs to be extended someday to
        // support arrow-nulls, we can work on that.
        if (child->buffers[0] != nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_column_by_index: validity buffer "
                "unsupported here");
        }

        const char* data = (char*)child->buffers[2];

        if (data == nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_column_by_index: null data buffer");
        }

        if (strcmp(arrow_schema->children[column_index]->format, "U") != 0) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_column_by_index: expected Arrow "
                "large_string");
        }
        uint64_t* offsets = (uint64_t*)child->buffers[1];

        int num_cells = (int)child->length;
        std::vector<std::string> retval(num_cells);
        for (int j = 0; j < num_cells; j++) {
            std::string e(&data[offsets[j]], &data[offsets[j + 1]]);
            retval[j] = e;
        }

        return retval;
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

    // Throws if the array and the schema don't have the same
    // recursive child-counts.
    static void _check_shapes(
        ArrowArray* arrow_array, ArrowSchema* arrow_schema);

    // Throws if the table doesn't have the column name.
    static int64_t _get_column_index_from_name(
        const ArrowTable& arrow_table, std::string column_name);

    static ArrowArray* _get_and_check_column(
        const ArrowTable& arrow_table,
        int64_t column_index,
        int64_t expected_n_buffers);

};  // class ArrowAdapter

};  // namespace tiledbsoma
#endif
