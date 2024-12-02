#ifndef ARROW_ADAPTER_H
#define ARROW_ADAPTER_H

#include <any>
#include <format>
#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

// https://arrow.apache.org/docs/format/CDataInterface.html
// https://arrow.apache.org/docs/format/Columnar.html#buffer-listing-for-each-layout
// https://arrow.apache.org/docs/format/CDataInterface.html#exporting-a-simple-int32-array

// A general developer note: in several places we have
//
//     template <typename T>
//     static sometype foo(T arg) {
//         if (std::is_same_v<T, std::string>) {
//             throw std::runtime_error(...);
//         }
//         }D...
//     }
//
//     static sometype foo_string(std::string arg) { ... }
//
// -- with explicit `_string` suffix -- rather than
//
//     template <typename T>
//     static sometype foo(T arg) ...
//
//     template <>
//     static sometype foo(std::string arg) ...
//
// We're aware of the former but we've found it a bit fiddly across systems and
// compiler versions -- namely, with the latter we find it tricky to always
// avoid the <typename T> variant being templated with std::string. It's simple,
// explicit, and robust to go the `_string` suffix route, and it's friendlier to
// current and future maintainers of this code.

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
     * @brief Create a an ArrowSchema from TileDB Dimension
     *
     * @return ArrowSchema
     */
    static std::unique_ptr<ArrowSchema> arrow_schema_from_tiledb_dimension(
        const Dimension& dimension);

    /**
     * @brief Create a an ArrowSchema from TileDB Attribute
     *
     * @return ArrowSchema
     */
    static std::unique_ptr<ArrowSchema> arrow_schema_from_tiledb_attribute(
        Attribute& attribute, const Context& ctx, const Array& tiledb_array);

    /**
     * @brief Get members of the TileDB Schema in the form of a PlatformConfig
     *
     * @return PlatformConfig
     */
    static PlatformConfig platform_config_from_tiledb_schema(
        ArraySchema tiledb_schema);

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
        const std::unique_ptr<ArrowSchema>& arrow_schema,
        const ArrowTable& index_column_info,
        std::string soma_type,
        bool is_sparse = true,
        PlatformConfig platform_config = PlatformConfig(),
        const ArrowTable& spatial_column_info = {
            std::unique_ptr<ArrowArray>(nullptr),
            std::unique_ptr<ArrowSchema>(nullptr)});

    /**
     * @brief Get a TileDB dimension from an Arrow schema.
     *
     * @return std::pair<Dimension, bool> The TileDB dimension.
     */
    static Dimension tiledb_dimension_from_arrow_schema(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        ArrowArray* array,
        std::string soma_type,
        std::string_view type_metadata,
        std::string prefix = std::string(),
        std::string suffix = std::string(),
        PlatformConfig platform_config = PlatformConfig());

    /**
     * @brief Get a TileDB attribute with its enumeration from an Arrow schema.
     *
     * @return std::pair<Attribute, std::optional<Enumeration>>
     */
    static std::pair<Attribute, std::optional<Enumeration>>
    tiledb_attribute_from_arrow_schema(
        std::shared_ptr<Context> ctx,
        ArrowSchema* arrow_schema,
        std::string_view type_metadata,
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
     * @brief Keystroke saver to determine whether Arrow type is of string,
     * large string, binary, or large binary type.
     *
     * @param const char* Arrow data format
     * @return bool Whether the Arrow type represents a string type
     */
    static bool arrow_is_var_length_type(const char* format);

    /**
     * @brief Get TileDB datatype from Arrow format string.
     *
     * @param datatype TileDB datatype.
     * @param arrow_dtype_metadata Additional datatype info. Useful for
     * differentiating between BLOB and WKB.
     * @return std::string_view Arrow format string.
     */
    static tiledb_datatype_t to_tiledb_format(
        std::string_view arrow_dtype,
        std::string_view arrow_dtype_metadata = {});

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
     * @brief Creates a nanoarrow ArrowSchema which accommodates
     * a varying number of columns.
     *
     * Note that the parents and children in nanoarrow are both of type
     * ArrowSchema. This constructs the parent and not the children.
     */
    static std::unique_ptr<ArrowSchema> make_arrow_schema_parent(
        int num_columns);

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

    static void log_make_arrow_array_child(ArrowArray* child);

    template <typename T>
    static ArrowArray* make_arrow_array_child_var(
        const std::pair<std::vector<T>, std::vector<T>>& pair) {
        std::vector<T> v = pair.first;
        v.insert(v.end(), pair.second.begin(), pair.second.end());
        ArrowArray* child = make_arrow_array_child<T>(v);
        return child;
    }

    static ArrowArray* make_arrow_array_child_string(
        const std::pair<std::string, std::string>& pair) {
        std::vector<std::string> v({pair.first, pair.second});
        return make_arrow_array_child_string(v);
    }

    static ArrowArray* make_arrow_array_child_binary() {
        // Use malloc here, not new, to match ArrowAdapter::release_array
        auto arrow_array = (ArrowArray*)malloc(sizeof(ArrowArray));

        ArrowArrayInitFromType(
            arrow_array, ArrowType::NANOARROW_TYPE_LARGE_BINARY);

        return arrow_array;
    }

    template <typename T>
    static ArrowArray* make_arrow_array_child(const std::vector<T>& v) {
        // We're aware of template-specialization wherein we can
        // make a separate `make_arrow_array_child` with `template <>` (not
        // `make_arrow_array_child_string`, as we've done). However,
        // we find it a bit fiddly across different compilers to force
        // the compiler to find the string variant. It's far more intuitive
        // for the non-expert developers (and maybe even for the experts),
        // and far more robust for any future maintainers, to explicitly
        // (a) have a separate `..._string` variant; (b) throw here
        // if callsites don't use t.
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

    template <typename T, size_t D>
    static std::vector<std::array<T, D>> get_table_column_by_name(
        const ArrowTable& arrow_table, std::string column_name) {
        std::vector<std::any> data = get_table_any_column_by_name<D>(
            arrow_table, column_name);
        std::vector<std::array<T, D>> result;

        for (size_t i = 0; i < data.size(); ++i) {
            result.push_back(std::any_cast<std::array<T, D>>(data[i]));
        }

        return result;
    }

    /**
     * Return a copy of the data in a specified column of an arrow table.
     * Complex column types are supported. The for each sub column are an
     * std::array<T, 2> casted as an std::any object.
     */
    template <size_t D>
    static std::vector<std::any> get_table_any_column_by_name(
        const ArrowTable& arrow_table, std::string column_name) {
        int64_t index = _get_column_index_from_name(arrow_table, column_name);
        return get_table_any_column_by_index<D>(arrow_table, index);
    }

    template <size_t D>
    static std::vector<std::any> get_table_any_column_by_index(
        const ArrowTable& arrow_table, int64_t column_index) {
        ArrowArray* arrow_array = arrow_table.first.get();
        ArrowSchema* arrow_schema = arrow_table.second.get();
        _check_shapes(arrow_array, arrow_schema);

        if (arrow_array->n_children == 0) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_any_column_by_index: expected "
                "non-leaf "
                "node");
        }

        if (arrow_schema->n_children <= column_index) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_any_column_by_index: Column index out "
                "fo bounds.");
        }

        std::vector<std::any> result;

        ArrowArray* selected_array = arrow_array->children[column_index];
        ArrowSchema* selected_schema = arrow_schema->children[column_index];

        // Complex domain
        if (selected_array->n_children != 0) {
            for (int64_t i = 0; i < selected_schema->n_children; ++i) {
                ArrowArray* array = selected_array->children[i];
                ArrowSchema* schema = selected_schema->children[i];

                result.push_back(get_table_any_column<D>(array, schema));
            }
        } else {
            result.push_back(
                get_table_any_column<D>(selected_array, selected_schema));
        }

        return result;
    }

    template <size_t D>
    static std::any get_table_any_column(
        ArrowArray* array, ArrowSchema* schema) {
        auto tdb_type = to_tiledb_format(schema->format, "");

        if (array->n_children != 0) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_any_column: expected leaf "
                "node");
        }

        if (array->length != D) {
            throw std::runtime_error(std::format(
                "ArrowAdapter::get_table_any_column: expected {} elements", D));
        }

        if (strcmp(schema->format, "u") == 0 ||
            strcmp(schema->format, "z") == 0 ||
            strcmp(schema->format, "U") == 0 ||
            strcmp(schema->format, "Z") == 0) {
            if (array->n_buffers != 3) {
                throw std::runtime_error(
                    "ArrowAdapter::get_table_any_column: expected three "
                    "buffers");
            }

            if (array->buffers[0] != nullptr) {
                throw std::runtime_error(
                    "ArrowAdapter::get_table_any_column: validity buffer "
                    "unsupported here");
            }
            if (array->buffers[1] == nullptr) {
                throw std::runtime_error(
                    "ArrowAdapter::get_table_any_column: null "
                    "offsets buffer");
            }
            if (array->buffers[2] == nullptr) {
                throw std::runtime_error(
                    "ArrowAdapter::get_table_any_column: null data "
                    "buffer");
            }
        } else {
            if (array->n_buffers != 2) {
                throw std::runtime_error(
                    "ArrowAdapter::get_table_any_column: expected two "
                    "buffers");
            }

            if (array->buffers[0] != nullptr) {
                throw std::runtime_error(
                    "ArrowAdapter::get_table_any_column: validity buffer "
                    "unsupported here");
            }
            if (array->buffers[1] == nullptr) {
                throw std::runtime_error(
                    "ArrowAdapter::get_table_any_column: null data buffer");
            }
        }

        switch (tdb_type) {
            case TILEDB_BOOL: {
                return std::make_any<std::array<bool, D>>(std::to_array(
                    (bool(&)[D])(*((bool*)array->buffers[1]))));
            }
            case TILEDB_UINT8: {
                return std::make_any<std::array<uint8_t, D>>(std::to_array(
                    (uint8_t(&)[D])(*((uint8_t*)array->buffers[1]))));
            }
            case TILEDB_UINT16: {
                return std::make_any<std::array<uint16_t, D>>(std::to_array(
                    (uint16_t(&)[D])(*((uint16_t*)array->buffers[1]))));
            }
            case TILEDB_UINT32: {
                return std::make_any<std::array<uint32_t, D>>(std::to_array(
                    (uint32_t(&)[D])(*((uint32_t*)array->buffers[1]))));
            }
            case TILEDB_UINT64: {
                return std::make_any<std::array<uint64_t, D>>(std::to_array(
                    (uint64_t(&)[D])(*((uint64_t*)array->buffers[1]))));
            }
            case TILEDB_INT8: {
                return std::make_any<std::array<int8_t, D>>(std::to_array(
                    (int8_t(&)[D])(*((int8_t*)array->buffers[1]))));
            }
            case TILEDB_INT16: {
                return std::make_any<std::array<int16_t, D>>(std::to_array(
                    (int16_t(&)[D])(*((int16_t*)array->buffers[1]))));
            }
            case TILEDB_INT32: {
                return std::make_any<std::array<int32_t, D>>(std::to_array(
                    (int32_t(&)[D])(*((int32_t*)array->buffers[1]))));
            }
            case TILEDB_DATETIME_YEAR:
            case TILEDB_DATETIME_MONTH:
            case TILEDB_DATETIME_WEEK:
            case TILEDB_DATETIME_DAY:
            case TILEDB_DATETIME_HR:
            case TILEDB_DATETIME_MIN:
            case TILEDB_DATETIME_SEC:
            case TILEDB_DATETIME_MS:
            case TILEDB_DATETIME_US:
            case TILEDB_DATETIME_NS:
            case TILEDB_DATETIME_PS:
            case TILEDB_DATETIME_FS:
            case TILEDB_DATETIME_AS:
            case TILEDB_INT64: {
                return std::make_any<std::array<int64_t, D>>(std::to_array(
                    (int64_t(&)[D])(*((int64_t*)array->buffers[1]))));
            }
            case TILEDB_FLOAT32: {
                return std::make_any<std::array<float_t, D>>(std::to_array(
                    (float_t(&)[D])(*((float_t*)array->buffers[1]))));
            }
            case TILEDB_FLOAT64: {
                return std::make_any<std::array<double_t, D>>(std::to_array(
                    (double_t(&)[D])(*((double_t*)array->buffers[1]))));
            }
            case TILEDB_STRING_ASCII:
            case TILEDB_STRING_UTF8:
            case TILEDB_CHAR:
            case TILEDB_BLOB:
            case TILEDB_GEOM_WKT:
            case TILEDB_GEOM_WKB: {
                if (strcmp(schema->format, "u") == 0 ||
                    strcmp(schema->format, "z") == 0) {
                    uint32_t* offsets = (uint32_t*)array->buffers[1];
                    char* data = (char*)array->buffers[2];

                    std::array<std::string, D> result;
                    for (size_t i = 0; i < D; ++i) {
                        result[i] = std::string(
                            &data[offsets[i]], &data[offsets[i + 1]]);
                    }

                    return std::make_any<std::array<std::string, D>>(result);
                } else if (
                    strcmp(schema->format, "U") == 0 ||
                    strcmp(schema->format, "Z") == 0) {
                    uint64_t* offsets = (uint64_t*)array->buffers[1];
                    char* data = (char*)array->buffers[2];

                    std::array<std::string, D> result;
                    for (size_t i = 0; i < D; ++i) {
                        result[i] = std::string(
                            &data[offsets[i]], &data[offsets[i + 1]]);
                    }

                    return std::make_any<std::array<std::string, D>>(result);
                } else {
                    throw std::runtime_error(std::format(
                        "ArrowAdapter::get_table_any_column: Unknown "
                        "schema format {}",
                        schema->format));
                }
            } break;
            default:
                throw std::runtime_error(std::format(
                    "ArrowAdapter::get_table_any_column: Unknown "
                    "datatype {}",
                    tiledb::impl::type_to_str(tdb_type)));
                break;
        }
    }

    static void set_current_domain_slot(
        tiledb_datatype_t type,
        const void* buff,
        NDRectangle& ndrect,
        std::string name);

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

    static json _get_attrs_filter_list_json(const ArraySchema& tiledb_schema);

    static json _get_dims_list_json(const ArraySchema& tiledb_schema);

    static json _get_filter_list_json(FilterList filter_list);

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

    static void _set_spatial_dimensions(
        std::map<std::string, Dimension>& dims,
        const ArrowTable& spatial_column_info,
        std::string_view type_metadata,
        std::string soma_type,
        std::shared_ptr<Context> ctx,
        PlatformConfig platform_config);
};  // class ArrowAdapter
};  // namespace tiledbsoma
#endif
