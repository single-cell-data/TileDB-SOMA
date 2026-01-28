/**
 * @file   soma_array.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This declares the SOMAArray class.
 */

#ifndef SOMA_ARRAY
#define SOMA_ARRAY

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <future>
#include <optional>
#include <span>
#include <utility>

#include "../utils/common.h"
#include "common/logging/logger.h"
#include "enums.h"
#include "soma_object.h"

namespace tiledb {
class Array;
class ArraySchema;
class ArraySchemaEvolution;
class Enumeration;
struct CurrentDomain;
}  // namespace tiledb

struct ArrowArray;
struct ArrowSchema;

namespace tiledbsoma {

class CoordinateValueFilters;
class SOMAColumn;
class ManagedQuery;
struct PlatformSchemaConfig;

using StatusAndReason = std::pair<bool, std::string>;
using ArrowTable = std::pair<managed_unique_ptr<ArrowArray>, managed_unique_ptr<ArrowSchema>>;

class SOMAArray : public SOMAObject {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAArray object at the given URI.
     *
     * @param ctx SOMAContext
     * @param uri URI to create the SOMAArray
     * @param schema TileDB ArraySchema
     * @param soma_type SOMADataFrame, SOMADenseNDArray, or
     * SOMASparseNDArray
     */
    static void create(
        std::shared_ptr<SOMAContext> ctx,
        std::string_view uri,
        tiledb::ArraySchema schema,
        std::string_view soma_type,
        std::optional<std::string_view> soma_schema = std::nullopt,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open an array at the specified URI and return SOMAArray
     * object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param platform_config Config parameter dictionary
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::unique_ptr<SOMAArray> SOMAArray
     */
    static std::unique_ptr<SOMAArray> open(
        OpenMode mode,
        std::string_view uri,
        std::map<std::string, std::string> platform_config = {},
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open an array at the specified URI and return SOMAArray
     * object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx SOMAContext
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::unique_ptr<SOMAArray> SOMAArray
     */
    static std::unique_ptr<SOMAArray> open(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMAArray object
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param platform_config Config parameter dictionary
     * @param timestamp Timestamp
     */
    SOMAArray(
        OpenMode mode,
        std::string_view uri,
        std::map<std::string, std::string> platform_config,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Construct a new SOMAArray object
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx SOMAContext
     * @param timestamp Timestamp
     */
    SOMAArray(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Construct a new SOMAArray from a TileDB array
     *
     * @param ctx SOMAContext
     * @param arr TileDB array.
     * @param timestamp Timestamp range the array was opened at.
     */
    SOMAArray(
        std::shared_ptr<SOMAContext> ctx, std::shared_ptr<tiledb::Array> arr, std::optional<TimestampRange> timestamp);

    SOMAArray(const SOMAArray& other) = default;
    SOMAArray(SOMAArray&&) = default;

    SOMAArray() = delete;
    virtual ~SOMAArray() = default;

    /**
     * @brief Get URI of the SOMAArray.
     *
     * @return std::string URI
     */
    const std::string uri() const;

    /**
     * @brief Get context of the SOMAArray.
     *
     * @return SOMAContext
     */
    inline std::shared_ptr<SOMAContext> ctx() {
        return ctx_;
    }

    /**
     * Open the SOMAArray object.
     *
     * @param mode read or write
     * @param timestamp Timestamp
     */
    void open(OpenMode mode, std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Returns a shared pointer of the internal TileDB array.
     */
    inline std::shared_ptr<tiledb::Array> tiledb_array() {
        return arr_;
    }

    /**
     * Creates a new ManagedQuery for this array.
     *
     * @param name Name of the array.
     */
    ManagedQuery create_managed_query(std::string_view name = "unnamed") const;

    /**
     * Creates a new ManagedQuery for this array.
     *
     * @param ctx SOMA context to use for the query.
     * @param name Name of the array.
     */
    ManagedQuery create_managed_query(std::shared_ptr<SOMAContext> query_ctx, std::string_view name = "unnamed") const;

    /**
     * Creates a new CoordinateValueFilters for this array.
     */
    CoordinateValueFilters create_coordinate_value_filter() const;

    /**
     * Close the SOMAArray object.
     */
    void close();

    /**
     * Check if the SOMAArray is open.
     *
     * @return bool true if open
     */
    bool is_open() const;

    /**
     * Get whether the SOMAArray was open in read or write mode.
     *
     * @return OpenMode
     */
    inline OpenMode mode() const {
        return soma_mode_;
    }

    /**
     * @brief Get the number of dimensions.
     *
     * @return uint64_t Number of dimensions.
     */
    uint64_t ndim() const;

    /**
     * @brief Get the name of each dimension.
     *
     * @return std::vector<std::string> Name of each dimension.
     */
    std::vector<std::string> dimension_names() const;

    /**
     * @brief Get the index columns.
     *
     * @return A vector of the index columns in order.
     */
    std::vector<std::shared_ptr<SOMAColumn>> index_columns() const;

    /**
     * @brief Sees if the array has a dimension of the given name.
     *
     * @return bool
     */
    bool has_dimension_name(std::string_view name) const;

    /**
     * @brief Get the name of each attribute.
     *
     * @return std::vector<std::string> Name of each attribute.
     */
    std::vector<std::string> attribute_names() const;

    /**
     * @brief Consolidates and vacuums fragment metadata and commit files.
     *
     * @param modes List of modes to apply. By default, apply to fragment_meta
     * and commits
     */
    void consolidate_and_vacuum(std::vector<std::string> modes = {"fragment_meta", "commits"});

    /**
     * @brief Delete cells from the array.
     *
     * @param coord_filter Coordinate value filter defining the coordinates to delete.
     */
    void delete_cells(const CoordinateValueFilters& coord_filters);

    /**
     * @brief Delete cells from the array.
     *
     * @param coord_filter Coordinate value filter defining the coordinates to delete.
     * @param value_filter Additional value filter to constrain the delete by.
     */
    void delete_cells(const CoordinateValueFilters& coord_filters, const QueryCondition& value_filter);

    /**
     * @brief Get the TileDB ArraySchema. This should eventually
     * be removed in lieu of arrow_schema below.
     *
     * @return std::shared_ptr<ArraySchema> Schema
     */
    std::shared_ptr<tiledb::ArraySchema> tiledb_schema() const {
        return schema_;
    }

    /**
     * @brief Get the Arrow schema of the array.
     *
     * @return std::unique_ptr<ArrowSchema> Schema
     */
    managed_unique_ptr<ArrowSchema> arrow_schema(bool downcast_dict_of_large_var = false) const;

    /**
     * @brief Get the Arrow schema of a column in the array.
     *
     * @return ArrowSchema* Schema
     */
    ArrowSchema* arrow_schema_for_column(std::string column_name, bool downcast_dict_of_large_var = false) const;

    /**
     * @brief Get members of the schema (capacity, allows_duplicates,
     * tile_order, cell_order, offsets_filters, validity_filters, attr filters,
     * and dim filters) in the form of a PlatformSchemaConfig.
     *
     * @return PlatformSchemaConfig
     */
    PlatformSchemaConfig schema_config_options() const;

    /**
     * Set metadata key-value items to an open array. The array must
     * opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be added. UTF-8 encodings
     *     are acceptable.
     * @param value_type The datatype of the value.
     * @param value_num The value may consist of more than one items of the
     *     same datatype. This argument indicates the number of items in the
     *     value component of the metadata.
     * @param value The metadata value in binary form.
     * @param force A boolean toggle to suppress internal checks, defaults to
     *     false.
     *
     * @note The writes will take effect only upon closing the array.
     */
    void set_metadata(
        const std::string& key,
        tiledb_datatype_t value_type,
        uint32_t value_num,
        const void* value,
        bool force = false);

    /**
     * Delete a metadata key-value item from an open array. The array must
     * be opened in WRITE mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be deleted.
     *
     * @note The writes will take effect only upon closing the array.
     *
     * @note If the key does not exist, this will take no effect
     *     (i.e., the function will not error out).
     */
    void delete_metadata(const std::string& key, bool force = false);

    /**
     * @brief Given a key, get the associated value datatype, number of
     * values, and value in binary form. The array must be opened in READ
     * mode, otherwise the function will error out.
     *
     * The value may consist of more than one items of the same datatype.
     * Keys that do not exist in the metadata will be return NULL for the value.
     *
     * **Example:**
     * @code{.cpp}
     * // Open the array for reading
     * tiledbsoma::SOMAArray soma_array = SOMAArray::open(TILEDB_READ,
     "s3://bucket-name/group-name");
     * tiledbsoma::MetadataValue meta_val = soma_array->get_metadata("key");
     * std::string key = std::get<MetadataInfo::key>(meta_val);
     * tiledb_datatype_t dtype = std::get<MetadataInfo::dtype>(meta_val);
     * uint32_t num = std::get<MetadataInfo::num>(meta_val);
     * const void* value = *((const
     int32_t*)std::get<MetadataInfo::value>(meta_val));
     * @endcode
     *
     * @param key The key of the metadata item to be retrieved. UTF-8
     * encodings are acceptable.
     * @return MetadataValue (std::tuple<std::string, tiledb_datatype_t,
     * uint32_t, const void*>)
     */
    std::optional<MetadataValue> get_metadata(const std::string& key);

    /**
     * Get a mapping of all metadata keys with its associated value datatype,
     * number of values, and value in binary form.
     *
     * @return std::map<std::string, MetadataValue>
     */
    std::map<std::string, MetadataValue> get_metadata();

    /**
     * Check if the key exists in metadata from an open array. The array
     * must be opened in READ mode, otherwise the function will error out.
     *
     * @param key The key of the metadata item to be checked. UTF-8
     * encodings are acceptable.
     * @return true if the key exists, else false.
     */
    bool has_metadata(const std::string& key);

    /**
     * Return then number of metadata items in an open array. The array must
     * be opened in READ mode, otherwise the function will error out.
     */
    uint64_t metadata_num() const;

    /**
     * Return optional timestamp pair SOMAArray was opened with.
     */
    std::optional<TimestampRange> timestamp();

    /**
     * Retrieves the enumeration values from the array's TileDB schema,
     * for specified column names. Throws if any of the column names
     * are not present in the schema, or if any of them are present but
     * none are for a non-enumerated column.
     *
     * @tparam std::vector<std::string> column names
     * @return ArrowTable with as many columns as the size of column_names.
     */
    ArrowTable get_enumeration_values(std::vector<std::string> column_names);

    /**
     * Retrieves the enumeration values from the array's TileDB schema,
     * for the specified column name. Throws if the column name is not
     * present in the schema, or if it is present but is for a non-enumerated
     * column.
     *
     * @tparam std::string column name
     * @return ArrowTable with one column
     */
    std::pair<ArrowArray*, ArrowSchema*> get_enumeration_values_for_column(std::string column_name);

    /**
     * Adds new values to enumeration columns.
     *
     * If deduplicate is `false`, the provided values for columns must be new,
     * unique values.
     *
     * @param values A mapping of column names to Arrow tables of enumeration
     * values.
     * @param deduplicate If set to false, new and existing values must be
     * disjoint for each given column.
     */
    void extend_enumeration_values(
        std::map<std::string, std::pair<ArrowSchema*, ArrowArray*>> values, bool deduplicate);

    /**
     * Exposed for testing purposes within this library.
     * Not for use by Python/R.
     */
    tiledb::CurrentDomain get_current_domain_for_test() const;

    /**
     * @brief Returns true if the array has a non-empty current domain, else
     * false.  Note that at the core level it's "current domain" for all arrays;
     * at the SOMA-API level it's "upgraded_shape" for SOMASparseNDArray and
     * SOMADenseNDArray, and "upgraded_domain" for SOMADataFrame; here
     * we use the core language and defer to Python/R to conform to
     * SOMA-API syntax.
     */
    bool has_current_domain() const;

    /**
     * Returns the SOMA domain in its entirety, as an Arrow table for return to
     * Python/R.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     *   - soma maxdomain is core domain
     *   - core current domain is not accessed at the soma level
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    inline ArrowTable get_soma_domain() {
        if (has_current_domain()) {
            return _get_core_domainish(Domainish::kind_core_current_domain);
        } else {
            return _get_core_domainish(Domainish::kind_core_domain);
        }
    }

    /**
     * Returns the SOMA maxdomain in its entirety, as an Arrow table for return
     * to Python/R.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     *   - soma maxdomain is core domain
     *   - core current domain is not accessed at the soma level
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    ArrowTable get_soma_maxdomain() {
        return _get_core_domainish(Domainish::kind_core_domain);
    }

    /**
     * Returns the core non-empty domain in its entirety, as an Arrow
     * table for return to Python/R.
     */
    ArrowTable get_non_empty_domain() {
        return _get_core_domainish(Domainish::kind_non_empty_domain);
    }

    /**
     * Code-dedupe helper for core domain, core current domain, and core
     * non-empty domain.
     */
    ArrowTable _get_core_domainish(enum Domainish which_kind);

    /**
     * @brief Get the total number of unique cells in the array. Equivalent to
     * the COUNT aggregate.
     *
     * @return uint64_t Total number of unique cells
     */
    uint64_t nnz();

    /**
     * @brief Get the total number of cells in all fragments. Does not account
     * for duplicates or deletes, and will therefore be an upper bound on the
     * actual unique cell count. See `nnz`.
     */
    uint64_t fragment_cell_count();

    /**
     * @brief Get the current capacity of each dimension.
     *
     * This applies to arrays all of whose dims are of type int64_t: this
     * includes SOMASparseNDArray and SOMADenseNDArray, and default-indexed
     * SOMADataFrame.
     *
     * At the TileDB-SOMA level we call this "shape". At the TileDB Core
     * storage level this maps to "current domain".
     *
     * Further, we map this single n to the pair (0, n-1) since core permits a
     * doubly inclusive pair (lo, hi) on each dimension slot.
     *
     * @return A vector with length equal to the number of dimensions; each
     * value in the vector is the capacity of each dimension.
     */
    std::vector<int64_t> shape();

    /**
     * @brief Get the maximum resizable capacity of each dimension.
     *
     * This applies to arrays all of whose dims are of type int64_t: this
     * includes SOMASparseNDArray and SOMADenseNDArray, and default-indexed
     * SOMADataFrame.
     *
     * At the TileDB-SOMA level we call this "maxshape". At the TileDB Core
     * storage level this maps to "domain".
     *
     * Further, we map this single n to the pair (0, n-1) since core permits a
     * doubly inclusive pair (lo, hi) on each dimension slot.
     *
     * @return A vector with length equal to the number of dimensions; each
     * value in the vector is the maximum capacity of each dimension.
     */
    std::vector<int64_t> maxshape();

    /**
     * This wires up to Python/R to tell a user if they can call resize() on an
     * array without error. For single arrays, they could just call resize() and
     * take their chances -- but for experiment-level resize (e.g. append mode)
     * it's crucial that we provide a can-we-do-them-all pass through all arrays
     * in the experiment before attempting any of them.
     *
     * On failure, returns false and an error string suitable for showing
     * to the user; on success, returns true and the empty string.
     *
     * Failure reasons: the requested shape's dimension-count doesn't match the
     * arrays; the array doesn't have a shape set (they must call
     * upgrade_shape), or the requested shape doesn't fit within the array's
     * existing core domain.
     */
    StatusAndReason can_resize(const std::vector<int64_t>& newshape, std::string function_name_for_messages) {
        return _can_set_shape_helper(newshape, true, function_name_for_messages);
    }

    /**
     * This wires up to Python/R to tell a user if they can call
     * upgrade_shape() on an array without error. For single dataframes,
     * they could just call upgrade_shape() and take their chances -- but for
     * experiment-level resize (e.g. append mode) it's crucial that we provide a
     * can-we-do-them-all pass through all arrays in the experiment before
     * attempting any of them.
     *
     * On failure, returns false and an error string suitable for showing
     * to the user; on success, returns true and the empty string.
     *
     * Failure reasons: the requested shape's dimension-count doesn't match the
     * arrays; the array already has a shape set (they must call resize), the
     * requested shape doesn't fit within the array's existing core domain, or
     * the requested shape is a downsize of the array's existing core current
     * domain.
     */
    StatusAndReason can_upgrade_shape(const std::vector<int64_t>& newshape, std::string function_name_for_messages) {
        return _can_set_shape_helper(newshape, false, function_name_for_messages);
    }

    /**
     * This is similar to can_upgrade_shape, but it's a can-we call
     * for resize_soma_joinid_shape.
     */
    StatusAndReason can_resize_soma_joinid_shape(int64_t newshape, std::string function_name_for_messages) {
        return _can_set_soma_joinid_shape_helper(newshape, true, function_name_for_messages);
    }

    /**
     * This is similar to can_upgrade_shape, but it's a can-we call
     * for upgrade_soma_joinid_shape.
     */
    StatusAndReason can_upgrade_soma_joinid_shape(int64_t newshape, std::string function_name_for_messages) {
        return _can_set_soma_joinid_shape_helper(newshape, false, function_name_for_messages);
    }

    /**
     * This is for SOMADataFrame.
     */
    StatusAndReason can_upgrade_domain(const ArrowTable& newdomain, std::string function_name_for_messages) {
        return _can_set_domain_helper(newdomain, false, function_name_for_messages);
    }

    /**
     * This is for SOMADataFrame.
     */
    StatusAndReason can_change_domain(const ArrowTable& newdomain, std::string function_name_for_messages) {
        return _can_set_domain_helper(newdomain, true, function_name_for_messages);
    }

    /**
     * @brief Resize the shape (what core calls "current domain") up to the
     * maxshape (what core calls "domain").
     *
     * This applies to arrays all of whose dims are of type int64_t: this
     * includes SOMASparseNDArray and SOMADenseNDArray, and default-indexed
     * SOMADataFrame.
     *
     * @return Nothing. Raises an exception if the resize would be a downsize,
     * which is not supported.
     */
    void resize(const std::vector<int64_t>& newshape, std::string function_name_for_messages) {
        _set_shape_helper(newshape, true, function_name_for_messages);
    }

    /**
     * @brief Given an old-style array without current domain, sets its
     * current domain. This is applicable only to arrays having all dims
     * of int64 type. Namely, all SparseNDArray/DenseNDArray, and
     * default-indexed DataFrame.
     */
    void upgrade_shape(const std::vector<int64_t>& newshape, std::string function_name_for_messages) {
        _set_shape_helper(newshape, false, function_name_for_messages);
    }

    /**
     * @brief Increases the tiledbsoma shape up to at most the maxshape,
     * resizing the soma_joinid dimension if it is a dimension.
     *
     * While SOMA SparseNDArray and DenseNDArray, along with default-indexed
     * DataFrame, have int64_t dims, non-default-indexed DataFrame objects need
     * not: it is only required that they have a dim _or_ an attr called
     * soma_joinid. If soma_joinid is one of the dims, it will be resized while
     * the others will be preserved. If soma_joinid is not one of the dims,
     * nothing will be changed, as nothing _needs_ to be changed.
     *
     * @return Throws if the requested shape exceeds the array's create-time
     * maxshape. Throws if the array does not have current-domain support.
     */
    void resize_soma_joinid_shape(int64_t newshape, std::string function_name_for_messages) {
        return _set_soma_joinid_shape_helper(newshape, true, function_name_for_messages);
    }

    /**
     * @brief Increases the tiledbsoma shape up to at most the maxshape,
     * resizing the soma_joinid dimension if it is a dimension.
     *
     * While SOMA SparseNDArray and DenseNDArray, along with default-indexed
     * DataFrame, have int64_t dims, non-default-indexed DataFrame objects need
     * not: it is only required that they have a dim _or_ an attr called
     * soma_joinid. If soma_joinid is one of the dims, it will be resized while
     * the others will be preserved. If soma_joinid is not one of the dims,
     * nothing will be changed, as nothing _needs_ to be changed.
     *
     * @return Throws if the requested shape exceeds the array's create-time
     * maxshape. Throws if the array does not have current-domain support.
     */
    void upgrade_soma_joinid_shape(int64_t newshape, std::string function_name_for_messages) {
        return _set_soma_joinid_shape_helper(newshape, false, function_name_for_messages);
    }

    /**
     * This is for SOMADataFrame.  While resize_soma_joinid_shape allows the
     * user to do up the soma_joinid domain slot, without needing to specify
     * the rest (which is the common operation for experiment-level resize)
     * this allows the full-generality resize-every-index-column case
     * (which only applies to variant-indexed/non-standard dataframes).
     */
    void change_domain(const ArrowTable& newdomain, std::string function_name_for_messages) {
        _set_domain_helper(newdomain, true, function_name_for_messages);
    }

    /**
     * This is for SOMADataFrame.  While upgrade_soma_joinid_shape allows the
     * user to do up the soma_joinid domain slot, without needing to specify
     * the rest (which is the common operation for experiment-level resize)
     * this allows the full-generality resize-every-index-column case
     * (which only applies to variant-indexed/non-standard dataframes).
     */
    void upgrade_domain(const ArrowTable& newdomain, std::string function_name_for_messages) {
        _set_domain_helper(newdomain, false, function_name_for_messages);
    }

    std::shared_ptr<SOMAColumn> get_column(std::string_view name) const;

    std::shared_ptr<SOMAColumn> get_column(std::size_t index) const;

   protected:
    static bool _exists(std::string_view uri, std::string_view soma_type, std::shared_ptr<SOMAContext> ctx);

    // These return the shape and maxshape slots for the soma_joinid dim, if
    // the array has one. These are important test-points and dev-internal
    // access-points, in particular, for the tiledbsoma-io experiment-level
    // resizer.
    std::optional<int64_t> _maybe_soma_joinid_shape();
    std::optional<int64_t> _maybe_soma_joinid_maxshape();

    /**
     * @brief Run a delete query on the array.
     *
     * @param delete_cond The query condition that specifies which cells to delete.
     */
    void delete_cells_impl(const QueryCondition& delete_cond);

    static tiledb::Array _create(
        std::shared_ptr<SOMAContext> ctx,
        std::string_view uri,
        ArraySchema schema,
        std::string_view soma_type,
        std::optional<std::string_view> soma_schema,
        std::optional<TimestampRange> timestamp);

    // SOMAArray URI
    std::string uri_;

    // NB: the Array dtor REQUIRES that this context be alive, so member
    // declaration order is significant.  Context (ctx_) MUST be declared
    // BEFORE Array (arr_) so that ctx_ will be destructed last.

    // SOMA context
    std::shared_ptr<SOMAContext> ctx_;

    // Array associated with SOMAArray
    std::shared_ptr<tiledb::Array> arr_;

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    /**
     * This can only be used when getting the enumeration for a column
     * which already exists, which is of enumerated type, and for which
     * the core enumeration already exists. Example uses including
     * getting the list of distinct enumeration values, or extending
     * the existing enumeration.
     */
    tiledb::Enumeration get_existing_enumeration_for_column(std::string column_name);

    /**
     * The caller must check the return value for .is_empty() to see if this is
     * a new-style array with current-domain support (.is_empty() is false) , or
     * an old-style array without current-domain support (.is_empty() is true).
     * We could implement this as a std::optional<CurrentDomain> return value
     * here, but, that would be a redundant indicator.
     */
    tiledb::CurrentDomain _get_current_domain() const;

    /**
     * This is a code-dedupe helper for can_resize and can_upgrade_shape.
     */
    StatusAndReason _can_set_shape_helper(
        const std::vector<int64_t>& newshape, bool must_already_have, std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for can_change_domain and can_upgrade_domain.
     */
    StatusAndReason _can_set_domain_helper(
        const ArrowTable& newdomain, bool must_already_have, std::string function_name_for_messages);

    /**
     * This is a second-level code-dedupe helper for _can_set_shape_helper.
     */
    StatusAndReason _can_set_shape_domainish_subhelper(
        const std::vector<int64_t>& newshape, std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper for can_upgrade_domain.
     */
    StatusAndReason _can_set_dataframe_domainish_subhelper(
        const ArrowTable& newdomain, std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper for can_resize_soma_joinid_shape and
     * can_upgrade_domain_soma_joinid_shape.
     */
    StatusAndReason _can_set_soma_joinid_shape_helper(
        int64_t newshape, bool must_already_have, std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for resize and upgrade_shape.
     */
    void _set_shape_helper(
        const std::vector<int64_t>& newshape, bool must_already_have, std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for resize_soma_joinid_shape and
     * upgrade_soma_joinid_shape.
     */
    void _set_soma_joinid_shape_helper(
        int64_t newshape, bool must_already_have, std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for change_domain and upgrade_domain.
     */
    void _set_domain_helper(
        const ArrowTable& newdomain, bool must_already_have, std::string function_name_for_messages);

    /**
     * While SparseNDArray, DenseNDArray, and default-indexed DataFrame
     * have int64 dims, variant-indexed DataFrames do not. This helper
     * lets us pre-check any attempts to treat dims as if they were int64.
     */
    bool _dims_are_int64();

    /**
     * Same, but throws.
     */
    void _check_dims_are_int64();

    /**
     * With old shape: core domain used to map to tiledbsoma shape; core current
     * domain did not exist.
     *
     * With new shape: core domain maps to tiledbsoma maxshape;
     * core current_domain maps to tiledbsoma shape.
     *
     * Here we distinguish between user-side API, and core-side implementation.
     */
    std::vector<int64_t> _shape_via_tiledb_domain();
    std::vector<int64_t> _shape_via_tiledb_current_domain();
    std::optional<int64_t> _maybe_soma_joinid_shape_via_tiledb_current_domain();
    std::optional<int64_t> _maybe_soma_joinid_shape_via_tiledb_domain();

    /**
     * Convenience function for creating an ArraySchemaEvolution object
     * referencing this array's context pointer, along with its open-at
     * timestamp (if any).
     */
    tiledb::ArraySchemaEvolution _make_se();

    // Array associated with metadata_. If open mode is not "read", then a second
    // array is opened. Otherwise, this is a reference to `arr_`. This needs to be
    // kept open to make metadata value pointer in the metadata cache accessible.
    std::shared_ptr<tiledb::Array> meta_cache_arr_;

    // Metadata cache
    std::map<std::string, MetadataValue> metadata_;

    // SOMAColumn list
    std::vector<std::shared_ptr<SOMAColumn>> columns_;

    // Read timestamp range (start, end)
    std::optional<TimestampRange> timestamp_;

    // The mode the SOMAArray is opened in.
    OpenMode soma_mode_;

    // The TileDB ArraySchema. The schema is inaccessible when the TileDB Array
    // is closed or opened in write mode which means we cannot use arr->schema()
    // directly in those cases. Here, we store a copy of the schema so that it
    // can be accessed in any mode
    std::shared_ptr<tiledb::ArraySchema> schema_;
};

}  // namespace tiledbsoma

#endif  // SOMA_ARRAY
