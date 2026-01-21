/**
 * @file   arrow_adapter.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines the ArrowAdapter class.
 */

#ifndef ARROW_ADAPTER_H
#define ARROW_ADAPTER_H

#include <any>
#include <concepts>
#include <future>
#include <span>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>
#include "../tiledb_adapter/platform_config.h"

// https://arrow.apache.org/docs/format/CDataInterface.html
// https://arrow.apache.org/docs/format/Columnar.html#buffer-listing-for-each-layout
// https://arrow.apache.org/docs/format/CDataInterface.html#exporting-a-simple-int32-array

#include "nanoarrow/nanoarrow.hpp"
#include "nlohmann/json.hpp"

namespace tiledbsoma {

using namespace tiledb;
using json = nlohmann::json;

class ReadColumnBuffer;
class ArrayBuffers;
class SOMACoordinateSpace;

/**
 * @brief The ArrowBuffer holds a shared pointer to a ColumnBuffer, which
 * manages the lifetime of a ColumnBuffer used to back an Arrow array.
 *
 * The ArrowArray.release callback will delete the ArrowBuffer, and
 * automatically decrement the use count of the ColumnBuffer's shared pointer.
 *
 */
// struct ArrowBuffer {
//     ArrowBuffer(CArrayColumnBuffer* buffer, bool large_offsets = true);
//     ArrowBuffer(ReadColumnBuffer* buffer, bool large_offsets = true);
//     ArrowBuffer(const Enumeration& enumeration, bool large_offsets = true);

//     std::span<std::byte> data_;
//     std::span<int64_t> large_offsets_;
//     std::span<int32_t> small_offsets_;
//     std::span<uint8_t> validity_;

//     size_t length;
//     std::string name;

//    private:
//     std::unique_ptr<std::byte[]> data_buffer_;
//     std::unique_ptr<std::byte[]> offset_buffer_;
//     std::unique_ptr<std::byte[]> validity_buffer_;
// };

// struct PrivateArrowBuffer {
//     PrivateArrowBuffer(const std::shared_ptr<ArrowBuffer>& buffer)
//         : buffer_(buffer) {
//     }

//     std::shared_ptr<ArrowBuffer> buffer_;
// };

template <typename T>
using managed_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

template <typename T, typename... Args>
    requires std::same_as<T, ArrowArray> || std::same_as<T, ArrowSchema>
managed_unique_ptr<T> make_managed_unique(Args&&... args) {
    return managed_unique_ptr<T>(new T(std::forward<Args>(args)...), [](T* arrow_struct) {
        if (arrow_struct->release != nullptr) {
            arrow_struct->release(arrow_struct);
        }

        delete arrow_struct;
    });
}

using ArrowTable = std::pair<managed_unique_ptr<ArrowArray>, managed_unique_ptr<ArrowSchema>>;

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
     * @brief Convert ArrayBuffer to an list of Arrow arrays.
     *
     * @return std::vector<std::pair<std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>>>
     */
    static std::vector<std::pair<managed_unique_ptr<ArrowArray>, managed_unique_ptr<ArrowSchema>>> buffer_to_arrow(
        std::shared_ptr<ArrayBuffers> buffers, bool downcast_dict_of_large_var = false);

    // /**
    //  * @brief Convert ColumnBuffer to an Arrow array.
    //  *
    //  * @return std::pair<std::unique_ptr<ArrowArray>,
    //  * std::unique_ptr<ArrowSchema>>
    //  */
    // static std::pair<managed_unique_ptr<ArrowArray>, managed_unique_ptr<ArrowSchema>> to_arrow(
    //     std::shared_ptr<ReadColumnBuffer> column,
    //     const std::unordered_map<std::string, std::shared_future<std::shared_ptr<ArrowBuffer>>>& enumerations,
    //     bool downcast_dict_of_large_var = false);

    /** @brief Create a an ArrowSchema from TileDB Dimension
     *
     * @return ArrowSchema
     */
    static ArrowSchema* arrow_schema_from_tiledb_dimension(const Dimension& dimension);

    /**
     * @brief Create a an ArrowSchema from TileDB Attribute
     *
     * @return ArrowSchema
     */
    static ArrowSchema* arrow_schema_from_tiledb_attribute(
        const Attribute& attribute, const Context& ctx, const Array& tiledb_array, bool downcast_dict_of_large_var);

    /**
     * @brief Create a TileDB ArraySchema from ArrowSchema and additional JSON
     * encoded metadata to handle higher level SOMA construct not supported by
     * TileDB (e.g. SOMAColumn)
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
     * @return std::tuple<tiledb::ArraySchema, nlohmann::json>
     */
    static std::tuple<ArraySchema, nlohmann::json> tiledb_schema_from_arrow_schema(
        std::shared_ptr<Context> ctx,
        const managed_unique_ptr<ArrowSchema>& arrow_schema,
        const ArrowTable& index_column_info,
        const std::optional<SOMACoordinateSpace>& coordinate_space,
        std::string soma_type,
        bool is_sparse = true,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<std::pair<int64_t, int64_t>> timestamp_range = std::nullopt);

    /**
     * @brief Get a TileDB attribute with its enumeration from an Arrow schema.
     *
     * @return std::pair<Attribute, std::optional<Enumeration>>
     */
    static std::pair<Attribute, std::optional<Enumeration>> tiledb_attribute_from_arrow_schema(
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
    static std::string_view to_arrow_format(tiledb_datatype_t tiledb_dtype, bool use_large = true);

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
    static tiledb_datatype_t to_tiledb_format(std::string_view arrow_dtype, std::string_view arrow_dtype_metadata = {});

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
    static managed_unique_ptr<ArrowSchema> make_arrow_schema(
        const std::vector<std::string>& names, const std::vector<tiledb_datatype_t>& tiledb_datatypes);

    /**
     * @brief Creates a nanoarrow ArrowSchema given a names and a TileDB
     * datatype.
     *
     * This constructs the child element, for a single column/attribute.
     */
    static ArrowSchema* make_arrow_schema_child(std::string name, tiledb_datatype_t tiledb_datatype);

    /**
     * @brief Creates a nanoarrow ArrowSchema which accommodates
     * a varying number of columns.
     *
     * Note that the parents and children in nanoarrow are both of type
     * ArrowSchema. This constructs the parent and not the children.
     */
    static managed_unique_ptr<ArrowSchema> make_arrow_schema_parent(
        size_t num_columns, std::string_view name = "parent");

    /**
     * @brief Creates a nanoarrow ArrowArray which accommodates
     * a varying number of columns.
     *
     * Note that the parents and children in nanoarrow are both of type
     * ArrowArray. This constructs the parent and not the children.
     */
    static managed_unique_ptr<ArrowArray> make_arrow_array_parent(size_t num_columns);

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
    static ArrowArray* make_arrow_array_child_var(const std::pair<std::vector<T>, std::vector<T>>& pair) {
        std::vector<T> v = pair.first;
        v.insert(v.end(), pair.second.begin(), pair.second.end());
        ArrowArray* child = make_arrow_array_child<T>(v);
        return child;
    }

    static ArrowArray* make_arrow_array_child_string(const std::pair<std::string, std::string>& pair) {
        std::vector<std::string> v({pair.first, pair.second});
        return make_arrow_array_child_string(v);
    }

    static ArrowArray* make_arrow_array_child_binary() {
        // Use malloc here, not new, to match ArrowAdapter::release_array
        auto arrow_array = (ArrowArray*)malloc(sizeof(ArrowArray));

        ArrowArrayInitFromType(arrow_array, ArrowType::NANOARROW_TYPE_LARGE_BINARY);

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
        // * Slot 1 is data, void* but will be dereferenced as T*
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
    static ArrowArray* make_arrow_array_child_string(const std::vector<std::string>& v) {
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
        // * Slot 2 is data, void* but will be dereferenced as T*
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

    // Given TileDB 8-bit booleans, packs them to 1-bit Arrow-style booleans.
    static ArrowArray* make_arrow_array_child_bool(const std::vector<uint8_t>& v) {
        // Use malloc here, not new, to match ArrowAdapter::release_array
        auto arrow_array = (ArrowArray*)malloc(sizeof(ArrowArray));

        size_t n = v.size();

        arrow_array->length = n;  // Number of bools, not number of bytes
        arrow_array->null_count = 0;
        arrow_array->offset = 0;

        // Two-buffer model for non-string data:
        // * Slot 0 is the Arrow validity buffer which we leave null
        // * Slot 1 is data, void* but will be dereferenced as T*
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
        // Use malloc here, not new, to match ArrowAdapter::release_array.
        // When computing the number of bits needed to store the n bits, round
        // up.
        auto nbytes = (n + 7) / 8;
        uint8_t* dest = (uint8_t*)malloc(nbytes * sizeof(uint8_t));
        memset(dest, 0, nbytes);

        for (size_t i = 0; i < n; i++) {
            if (v[i]) {
                ArrowBitSet(dest, i);
            }
        }

        arrow_array->buffers[1] = (void*)dest;

        log_make_arrow_array_child(arrow_array);

        return arrow_array;
    }

    // These table-column getters are, as of this writing, intended
    // primarily for keystroke-reduction in unit-test cases.

    template <typename T>
    static std::vector<T> get_table_non_string_column_by_name(const ArrowTable& arrow_table, std::string column_name) {
        int64_t index = _get_column_index_from_name(arrow_table, column_name);
        return get_table_non_string_column_by_index<T>(arrow_table, index);
    }

    static std::vector<std::string> get_table_string_column_by_name(
        const ArrowTable& arrow_table, std::string column_name) {
        int64_t index = _get_column_index_from_name(arrow_table, column_name);
        return get_table_string_column_by_index(arrow_table, index);
    }

    /**
     * Returns a copy of the data in a specified non-string column of an
     * ArrowTable as a standard/non-Arrow C++ object.
     */
    template <typename T>
    static std::vector<T> get_table_non_string_column_by_index(const ArrowTable& arrow_table, int64_t column_index) {
        ArrowArray* arrow_array = arrow_table.first.get();
        ArrowSchema* arrow_schema = arrow_table.second.get();
        _check_shapes(arrow_array, arrow_schema);

        if (std::is_same_v<T, std::string>) {
            throw std::runtime_error(
                "SOMAArray::get_table_non_string_column_by_index: "
                "template-specialization failure.");
        }

        ArrowArray* child_array = _get_and_check_column(arrow_table, column_index, 2);

        return get_array_non_string_column<T>(child_array);
    }

    /**
     * Returns a copy of the data in a specified string column of an
     * ArrowTable as a standard/non-Arrow C++ object.
     */
    static std::vector<std::string> get_table_string_column_by_index(
        const ArrowTable& arrow_table, int64_t column_index) {
        ArrowArray* arrow_array = arrow_table.first.get();
        ArrowSchema* arrow_schema = arrow_table.second.get();
        _check_shapes(arrow_array, arrow_schema);

        ArrowArray* child_array = _get_and_check_column(arrow_table, column_index, 3);

        const ArrowSchema* child_schema = arrow_schema->children[column_index];

        return get_array_string_column(child_array, child_schema);
    }

    /**
     * Returns a copy of the data in a specified non-string column of an
     * ArrowArray as a standard/non-Arrow C++ object. There is no
     * ArrowSchema* argument, as the caller must have determined the Arrow
     * type, inferred a C++ type, and have invoked this method with the
     * appropriate C++ type.
     *
     * This is a helper for get_table_non_string_column_by_index; also
     * exposed for callsites which have access to child objects which are
     * not top-level ArrowTables.
     */
    template <typename T>
    static std::vector<T> get_array_non_string_column(const ArrowArray* arrow_array) {
        if (arrow_array->n_children != 0) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_non_string_column: expected leaf "
                "node");
        }
        if (arrow_array->n_buffers != 2) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_non_string_column: expected two "
                "buffers");
        }

        if (std::is_same_v<T, std::string>) {
            throw std::runtime_error(
                "SOMAArray::get_array_non_string_column: "
                "template-specialization "
                "failure.");
        }

        // Two-buffer model for non-string data:
        // * Slot 0 is the Arrow validity buffer which we leave null
        // * Slot 1 is data, void* but will be dereferenced as T*
        // * There is no offset information

        // For our purposes -- reporting domains, etc. -- we don't use the
        // Arrow validity buffers. If this class needs to be extended
        // someday to support arrow-nulls, we can work on that.
        if (arrow_array->buffers[0] != nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_non_string_column: validity "
                "buffer "
                "unsupported here");
        }
        if (arrow_array->buffers[1] == nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_non_string_column: null data "
                "buffer");
        }

        const void* vdata = arrow_array->buffers[1];
        if (vdata == nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_non_string_column: null data "
                "buffer");
        }

        const T* data = (T*)vdata;

        std::vector<T> retval(arrow_array->length);
        for (auto i = 0; i < arrow_array->length; i++) {
            retval[i] = data[i];
        }
        return retval;
    }

    /**
     * Returns a copy of the data in a specified non-string column of an
     * ArrowArray as a standard/non-Arrow C++ object.
     * Helper for get_table_string_column_by_index. Also exposed for
     * callsites which have access to child objects which are not top-level
     * ArrowTables.
     */
    static std::vector<std::string> get_array_string_column(
        const ArrowArray* arrow_array, const ArrowSchema* arrow_schema) {
        if (arrow_array->n_children != 0 || arrow_schema->n_children != 0) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_string_column: expected leaf "
                "node");
        }
        if (arrow_array->n_buffers != 3) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_string_column: expected three "
                "buffers");
        }

        // Three-buffer model for string data:
        // * Slot 0 is the Arrow uint8_t* validity buffer
        // * Slot 1 is the Arrow offsets buffer: uint32_t* for Arrow string
        //   or uint64_t* for Arrow large_string
        // * Slot 2 is data, void* but will be dereferenced as T*

        // For our purposes -- reporting domains, etc. -- we don't use the
        // Arrow validity buffers. If this class needs to be extended
        // someday to support arrow-nulls, we can work on that.
        if (arrow_array->buffers[0] != nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_string_column: validity buffer "
                "unsupported here");
        }
        if (arrow_array->buffers[1] == nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_string_column: null "
                "offsets buffer");
        }
        if (arrow_array->buffers[2] == nullptr) {
            throw std::runtime_error(
                "ArrowAdapter::get_array_string_column: null data "
                "buffer");
        }

        const char* data = (char*)arrow_array->buffers[2];

        if (strcmp(arrow_schema->format, "u") == 0 || strcmp(arrow_schema->format, "z") == 0) {
            uint32_t* offsets = (uint32_t*)arrow_array->buffers[1];
            int num_cells = (int)arrow_array->length;
            std::vector<std::string> retval(num_cells);
            for (int j = 0; j < num_cells; j++) {
                std::string e(&data[offsets[j]], &data[offsets[j + 1]]);
                retval[j] = e;
            }
            return retval;

        } else if (strcmp(arrow_schema->format, "U") == 0 || strcmp(arrow_schema->format, "Z") == 0) {
            uint64_t* offsets = (uint64_t*)arrow_array->buffers[1];
            int num_cells = (int)arrow_array->length;
            std::vector<std::string> retval(num_cells);
            for (int j = 0; j < num_cells; j++) {
                std::string e(&data[offsets[j]], &data[offsets[j + 1]]);
                retval[j] = e;
            }
            return retval;

        } else {
            throw std::runtime_error(
                "ArrowAdapter::get_array_string_column: expected "
                "Arrow string, large_string, binary, or large_binary");
        }
    }

    /**
     * Return a copy of the data in a specified column of an Arrow table.
     * Complex column types are supported. The type for each subcolumn is an
     * std::array<T, 2> casted as an std::any object.
     *
     * @tparam S The number of elements to retrieve
     */
    template <size_t S>
    static std::vector<std::any> get_table_any_column_by_name(
        const ArrowTable& arrow_table, std::string column_name, size_t offset) {
        int64_t index = _get_column_index_from_name(arrow_table, column_name);
        return get_table_any_column_by_index<S>(arrow_table, index, offset);
    }

    template <size_t S>
    static std::vector<std::any> get_table_any_column_by_index(
        const ArrowTable& arrow_table, int64_t column_index, size_t offset) {
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
                "ArrowAdapter::get_table_any_column_by_index: column index "
                "out "
                "of bounds.");
        }

        std::vector<std::any> result;

        ArrowArray* selected_array = arrow_array->children[column_index];
        ArrowSchema* selected_schema = arrow_schema->children[column_index];

        if (strcmp(selected_schema->format, "+s") == 0) {
            // Complex domain like the ones required by `GeometryColumn`
            // expect a struct containing lists of values
            for (int64_t i = 0; i < selected_schema->n_children; ++i) {
                ArrowArray* array = selected_array->children[i];
                ArrowSchema* schema = selected_schema->children[i];

                result.push_back(get_table_any_column<S>(array, schema, offset));
            }
        } else {
            result.push_back(get_table_any_column<S>(selected_array, selected_schema, offset));
        }

        return result;
    }

    /**
     * Read a part of an Arrow array to an std::array and cast into an
     * std::any object.
     *
     * @example get_table_any_column<3>(array, schema, offset) will return
     * an std::array<T, 3> where T is the appropriate type based on the
     * Arrow array format casted as an std::any object. The std::array will
     * skip as many elements as specified by `offset` and copy the next 3
     * from the Arrow array.
     *
     * @tparam S The number of elements to read.
     *
     * @param array The Arrow array to read the data from. The array should
     * be a leaf node with no subarrays.
     * @param schema The Arrow schema of the given array.
     * @param offset The number of elements to skip from the beginning of
     * the Arrow array
     *
     * @remarks This method's usage is to extract specific subranges of
     * ArrowArray data and they come in handy during ArrowSchema ->
     * TileDBSchema where the Arrow array provided has 5 values per
     * dimension and we only need the last 2 to set the current domain.
     *
     * `S` is required to be a template parameter to specify the
     * compile-time size of the underlying std::array that will hold the
     * extracted data.
     *
     * As to using std::variant, adding more SOMAColumn types would require
     * changing multiple variants. The use of std::any here is to enable
     * runtime polymorphism and indirectly introduces a runtime type check
     * (via any_cast, make_any) between the templated function and the
     * actual dimension type. std::variant can provide all the above; this
     * is a stylistic choice.
     */
    template <size_t S>
    static std::any get_table_any_column(ArrowArray* array, ArrowSchema* schema, size_t offset) {
        auto tdb_type = to_tiledb_format(schema->format, "");

        if (array->n_children != 0) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_any_column: expected leaf "
                "node");
        }

        if (array->length < static_cast<int64_t>(S + offset)) {
            throw std::runtime_error(
                "ArrowAdapter::get_table_any_column: expected at least " + std::to_string(S + offset) + " elements");
        }

        if (strcmp(schema->format, "u") == 0 || strcmp(schema->format, "z") == 0 || strcmp(schema->format, "U") == 0 ||
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
                throw std::runtime_error("ArrowAdapter::get_table_any_column: null data buffer");
            }
        }

        size_t arrow_offset = static_cast<size_t>(array->offset);

        switch (tdb_type) {
            case TILEDB_BOOL: {
                std::array<bool, S> result;
                for (size_t i = 0; i < S; ++i) {
                    result[i] = static_cast<bool>(
                        ArrowBitGet(static_cast<const uint8_t*>(array->buffers[1]), i + arrow_offset + offset));
                }
                return std::make_any<std::array<bool, S>>(result);
            }
            case TILEDB_UINT8:
                return std::make_any<std::array<uint8_t, S>>(
                    std::to_array((uint8_t (&)[S])(*((uint8_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_UINT16:
                return std::make_any<std::array<uint16_t, S>>(
                    std::to_array((uint16_t (&)[S])(*((uint16_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_UINT32:
                return std::make_any<std::array<uint32_t, S>>(
                    std::to_array((uint32_t (&)[S])(*((uint32_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_UINT64:
                return std::make_any<std::array<uint64_t, S>>(
                    std::to_array((uint64_t (&)[S])(*((uint64_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_INT8:
                return std::make_any<std::array<int8_t, S>>(
                    std::to_array((int8_t (&)[S])(*((int8_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_INT16:
                return std::make_any<std::array<int16_t, S>>(
                    std::to_array((int16_t (&)[S])(*((int16_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_INT32:
                return std::make_any<std::array<int32_t, S>>(
                    std::to_array((int32_t (&)[S])(*((int32_t*)array->buffers[1] + arrow_offset + offset))));
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
            case TILEDB_INT64:
                return std::make_any<std::array<int64_t, S>>(
                    std::to_array((int64_t (&)[S])(*((int64_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_FLOAT32:
                return std::make_any<std::array<float_t, S>>(
                    std::to_array((float_t(&)[S])(*((float_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_FLOAT64:
                return std::make_any<std::array<double_t, S>>(
                    std::to_array((double_t(&)[S])(*((double_t*)array->buffers[1] + arrow_offset + offset))));
            case TILEDB_STRING_ASCII:
            case TILEDB_STRING_UTF8:
            case TILEDB_CHAR:
            case TILEDB_GEOM_WKT: {
                if (strcmp(schema->format, "u") == 0 || strcmp(schema->format, "z") == 0) {
                    auto offsets = static_cast<const uint32_t*>(array->buffers[1]);
                    auto data = static_cast<const char*>(array->buffers[2]);

                    std::array<std::string, S> result;
                    for (size_t i = arrow_offset + offset; i < arrow_offset + S + offset; ++i) {
                        if (offsets[i + 1] - offsets[i] != 0) {
                            result[i - arrow_offset - offset] = std::string(&data[offsets[i]], &data[offsets[i + 1]]);
                        }
                    }

                    return std::make_any<std::array<std::string, S>>(result);
                } else if (strcmp(schema->format, "U") == 0 || strcmp(schema->format, "Z") == 0) {
                    auto offsets = static_cast<const uint64_t*>(array->buffers[1]);
                    auto data = static_cast<const char*>(array->buffers[2]);

                    std::array<std::string, S> result;
                    for (size_t i = arrow_offset + offset; i < arrow_offset + S + offset; ++i) {
                        if (offsets[i + 1] - offsets[i] != 0) {
                            result[i - arrow_offset - offset] = std::string(&data[offsets[i]], &data[offsets[i + 1]]);
                        }
                    }

                    return std::make_any<std::array<std::string, S>>(result);
                } else {
                    throw std::runtime_error(
                        "ArrowAdapter::get_table_any_column: Unknown "
                        "schema format '" +
                        std::string(schema->format) + "'");
                }
            } break;
            case TILEDB_BLOB:
            case TILEDB_GEOM_WKB: {
                if (strcmp(schema->format, "u") == 0 || strcmp(schema->format, "z") == 0) {
                    auto offsets = static_cast<const uint32_t*>(array->buffers[1]);
                    auto data = static_cast<const std::byte*>(array->buffers[2]);

                    std::array<std::vector<std::byte>, S> result;
                    for (size_t i = arrow_offset + offset; i < arrow_offset + S + offset; ++i) {
                        if (offsets[i + 1] - offsets[i] != 0) {
                            std::copy(
                                &data[offsets[i]], &data[offsets[i + 1]], result[i - arrow_offset - offset].begin());
                        }
                    }

                    return std::make_any<std::array<std::vector<std::byte>, S>>(result);
                } else if (strcmp(schema->format, "U") == 0 || strcmp(schema->format, "Z") == 0) {
                    auto offsets = static_cast<const uint64_t*>(array->buffers[1]);
                    auto data = static_cast<const std::byte*>(array->buffers[2]);

                    std::array<std::vector<std::byte>, S> result;
                    for (size_t i = arrow_offset + offset; i < arrow_offset + S + offset; ++i) {
                        if (offsets[i + 1] - offsets[i] != 0) {
                            std::copy(
                                &data[offsets[i]], &data[offsets[i + 1]], result[i - arrow_offset - offset].begin());
                        }
                    }

                    return std::make_any<std::array<std::vector<std::byte>, S>>(result);
                } else {
                    throw std::runtime_error(
                        "ArrowAdapter::get_table_any_column: Unknown "
                        "schema format '" +
                        std::string(schema->format) + "'");
                }
            } break;
            default:
                throw std::runtime_error(
                    "ArrowAdapter::get_table_any_column: Unknown "
                    "datatype '" +
                    tiledb::impl::type_to_str(tdb_type) + "'");
                break;
        }
    }

    static managed_unique_ptr<ArrowArray> arrow_array_insert_at_index(
        managed_unique_ptr<ArrowArray> parent_array,
        std::vector<managed_unique_ptr<ArrowArray>> child_arrays,
        int64_t index);

    static managed_unique_ptr<ArrowSchema> arrow_schema_insert_at_index(
        managed_unique_ptr<ArrowSchema> parent_schema,
        std::vector<managed_unique_ptr<ArrowSchema>> child_schemas,
        int64_t index);

    static managed_unique_ptr<ArrowArray> arrow_array_remove_at_index(
        managed_unique_ptr<ArrowArray> array, int64_t index);

    static managed_unique_ptr<ArrowSchema> arrow_schema_remove_at_index(
        managed_unique_ptr<ArrowSchema> schema, int64_t index);

   private:
    static size_t _set_var_dictionary_buffers(
        Enumeration& enumeration, const Context& ctx, const void** buffers, bool downcast_dict_of_large_var = false);

    static size_t _set_dictionary_buffers(Enumeration& enumeration, const Context& ctx, const void** buffers);

    static size_t _set_bool_dictionary_buffers(Enumeration& enumeration, const Context& ctx, const void** buffers);

    static Dimension _create_dim(
        tiledb_datatype_t type, std::string name, const void* buff, std::shared_ptr<Context> ctx);

    // Throws if the array and the schema don't have the same
    // recursive child-counts.
    static void _check_shapes(ArrowArray* arrow_array, ArrowSchema* arrow_schema);

    // Throws if the table doesn't have the column name.
    static int64_t _get_column_index_from_name(const ArrowTable& arrow_table, std::string column_name);

    static ArrowArray* _get_and_check_column(
        const ArrowTable& arrow_table, int64_t column_index, int64_t expected_n_buffers);

};  // class ArrowAdapter
};  // namespace tiledbsoma
#endif
