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

#include <format>
#include <future>
#include <span>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>
#include "../utils/arrow_adapter.h"
#include "enums.h"
#include "logger_public.h"
#include "managed_query.h"
#include "soma_column.h"
#include "soma_object.h"

// ================================================================
// Some general developer notes:
//
// ----------------------------------------------------------------
// In several places we have
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
//
// ----------------------------------------------------------------
// These are several methods here for use nominally by SOMADataFrame. These
// could be moved in their entirety to SOMADataFrame, but that would entail
// moving several SOMAArray attributes from private to protected, which has
// knock-on effects on the order of constructor initializers, etc.: in total
// it's simplest to place these here and have SOMADataFrame invoke them.
// ================================================================

namespace tiledbsoma {
using namespace tiledb;

using StatusAndReason = std::pair<bool, std::string>;

class SOMAArray : public SOMAObject {
   public:
    friend class ManagedQuery;

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
        ArraySchema schema,
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

    SOMAArray(const SOMAArray& other)
        // Ensure protected attributes initalized first in a consistent ordering
        : uri_(other.uri_)
        , ctx_(other.ctx_)
        , arr_(other.arr_)
        // Initialize private attributes next to control the order of
        // destruction
        , metadata_(other.metadata_)
        , timestamp_(other.timestamp_)
        , schema_(other.schema_)
        , meta_cache_arr_(other.meta_cache_arr_) {
        fill_metadata_cache(timestamp_);
        fill_columns();
    }

    SOMAArray(
        std::shared_ptr<SOMAContext> ctx,
        std::shared_ptr<Array> arr,
        std::optional<TimestampRange> timestamp);

    SOMAArray(SOMAArray&&) = default;

    SOMAArray(const SOMAObject& other)
        : SOMAObject(other) {
    }

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
    std::shared_ptr<SOMAContext> ctx();

    /**
     * Open the SOMAArray object.
     *
     * @param mode read or write
     * @param timestamp Timestamp
     */
    void open(
        OpenMode mode, std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Return a new SOMAArray with the given mode at the current Unix timestamp.
     *
     * @param mode if the OpenMode is not given, If the SOMAObject was opened in
     * READ mode, reopen it in WRITE mode and vice versa
     * @param timestamp Timestamp
     */
    std::unique_ptr<SOMAArray> reopen(
        OpenMode mode, std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * Close the SOMAArray object.
     */
    void close();

    /**
     * Check if the SOMAArray is open.
     *
     * @return bool true if open
     */
    bool is_open() const {
        return arr_->is_open();
    }

    /**
     * Get whether the SOMAArray was open in read or write mode.
     *
     * @return OpenMode
     */
    OpenMode mode() const {
        return arr_->query_type() == TILEDB_READ ? OpenMode::read :
                                                   OpenMode::write;
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
    void consolidate_and_vacuum(
        std::vector<std::string> modes = {"fragment_meta", "commits"});

    /**
     * @brief Get the TileDB ArraySchema. This should eventually
     * be removed in lieu of arrow_schema below.
     *
     * @return std::shared_ptr<ArraySchema> Schema
     */
    std::shared_ptr<ArraySchema> tiledb_schema() const {
        return schema_;
    }

    /**
     * @brief Get the Arrow schema of the array.
     *
     * @return std::unique_ptr<ArrowSchema> Schema
     */
    std::unique_ptr<ArrowSchema> arrow_schema() const {
        auto schema = ArrowAdapter::make_arrow_schema_parent(columns_.size());

        for (size_t i = 0; i < columns_.size(); ++i) {
            schema->children[i] = columns_[i]->arrow_schema_slot(*ctx_, *arr_);
        }

        return schema;
    }

    /**
     * @brief Get members of the schema (capacity, allows_duplicates,
     * tile_order, cell_order, offsets_filters, validity_filters, attr filters,
     * and dim filters) in the form of a PlatformSchemaConfig.
     *
     * @return PlatformSchemaConfig
     */
    PlatformSchemaConfig schema_config_options() const {
        return ArrowAdapter::platform_schema_config_from_tiledb(*schema_);
    }

    /**
     * @brief Get members of the schema (capacity, allows_duplicates,
     * tile_order, cell_order, offsets_filters, validity_filters, attr filters,
     * and dim filters) in the form of a PlatformConfig
     *
     * @return PlatformConfig
     */
    PlatformConfig config_options_from_schema() const {
        return ArrowAdapter::platform_config_from_tiledb_schema(*schema_);
    }

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
     * Validates input parameters before opening array.
     */
    void validate(OpenMode mode, std::optional<TimestampRange> timestamp);

    /**
     * Return optional timestamp pair SOMAArray was opened with.
     */
    std::optional<TimestampRange> timestamp();

    /**
     * Retrieves the non-empty domain from the array. This is the union of the
     * non-empty domains of the array fragments. Returns (0, 0) for empty
     * domains.
     */
    template <typename T>
    std::pair<T, T> non_empty_domain_slot(std::string_view name) const {
        return get_column(name)->non_empty_domain_slot<T>(*arr_);
    }

    /**
     * Retrieves the non-empty domain from the array. This is the union of the
     * non-empty domains of the array fragments. Return std::nullopt for empty
     * domains.
     */
    template <typename T>
    std::optional<std::pair<T, T>> non_empty_domain_slot_opt(
        std::string_view name) const {
        return get_column(name)->non_empty_domain_slot_opt<T>(*ctx_, *arr_);
    }

    /**
     * Exposed for testing purposes within this library.
     * Not for use by Python/R.
     */
    CurrentDomain get_current_domain_for_test() const {
        return _get_current_domain();
    }

    /**
     * @brief Returns true if the array has a non-empty current domain, else
     * false.  Note that at the core level it's "current domain" for all arrays;
     * at the SOMA-API level it's "upgraded_shape" for SOMASparseNDArray and
     * SOMADenseNDArray, and "upgraded_domain" for SOMADataFrame; here
     * we use the core language and defer to Python/R to conform to
     * SOMA-API syntax.
     */
    bool has_current_domain() const {
        return !_get_current_domain().is_empty();
    }

    /**
     * Returns the core current domain at the given dimension.
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
    template <typename T>
    std::pair<T, T> _core_current_domain_slot(std::string_view name) const {
        return get_column(name)->core_current_domain_slot<T>(*ctx_, *arr_);
    }

    /**
     * Returns the core domain at the given dimension.
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
    template <typename T>
    std::pair<T, T> _core_domain_slot(std::string_view name) const {
        return get_column(name)->core_domain_slot<T>();
    }

    /**
     * Returns the SOMA domain at the given dimension.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     */
    template <typename T>
    std::pair<T, T> soma_domain_slot(std::string_view name) const {
        if (has_current_domain()) {
            return _core_current_domain_slot<T>(name);
        } else {
            return _core_domain_slot<T>(name);
        }
    }

    /**
     * Returns the SOMA maxdomain at the given dimension.
     *
     * o For arrays with core current-domain support:
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma maxdomain is core domain
     */
    template <typename T>
    std::pair<T, T> soma_maxdomain_slot(std::string_view name) const {
        return _core_domain_slot<T>(name);
    }

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
    ArrowTable get_soma_domain() {
        if (has_current_domain()) {
            return _get_core_current_domain();
        } else {
            return _get_core_domain();
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
        return _get_core_domain();
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
     * This enables some code deduplication between core domain, core current
     * domain, and core non-empty domain.
     */
    template <typename T>
    std::pair<T, T> _core_domainish_slot(
        std::string_view name, enum Domainish which_kind) const {
        return get_column(name)->domain_slot<T>(*ctx_, *arr_, which_kind);
    }

    /**
     * @brief Get the total number of unique cells in the array.
     *
     * @return uint64_t Total number of unique cells
     */
    uint64_t nnz();

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
    StatusAndReason can_resize(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages) {
        return _can_set_shape_helper(
            newshape, true, function_name_for_messages);
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
    StatusAndReason can_upgrade_shape(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages) {
        return _can_set_shape_helper(
            newshape, false, function_name_for_messages);
    }

    /**
     * This is similar to can_upgrade_shape, but it's a can-we call
     * for resize_soma_joinid_shape.
     */
    StatusAndReason can_resize_soma_joinid_shape(
        int64_t newshape, std::string function_name_for_messages) {
        return _can_set_soma_joinid_shape_helper(
            newshape, true, function_name_for_messages);
    }

    /**
     * This is similar to can_upgrade_shape, but it's a can-we call
     * for upgrade_soma_joinid_shape.
     */
    StatusAndReason can_upgrade_soma_joinid_shape(
        int64_t newshape, std::string function_name_for_messages) {
        return _can_set_soma_joinid_shape_helper(
            newshape, false, function_name_for_messages);
    }

    /**
     * This is for SOMADataFrame.
     */
    StatusAndReason can_upgrade_domain(
        const ArrowTable& newdomain, std::string function_name_for_messages) {
        return _can_set_domain_helper(
            newdomain, false, function_name_for_messages);
    }

    /**
     * This is for SOMADataFrame.
     */
    StatusAndReason can_change_domain(
        const ArrowTable& newdomain, std::string function_name_for_messages) {
        return _can_set_domain_helper(
            newdomain, true, function_name_for_messages);
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
    void resize(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages) {
        _set_shape_helper(newshape, true, function_name_for_messages);
    }

    /**
     * @brief Given an old-style array without current domain, sets its
     * current domain. This is applicable only to arrays having all dims
     * of int64 type. Namely, all SparseNDArray/DenseNDArray, and
     * default-indexed DataFrame.
     */
    void upgrade_shape(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages) {
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
    void resize_soma_joinid_shape(
        int64_t newshape, std::string function_name_for_messages) {
        return _set_soma_joinid_shape_helper(
            newshape, true, function_name_for_messages);
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
    void upgrade_soma_joinid_shape(
        int64_t newshape, std::string function_name_for_messages) {
        return _set_soma_joinid_shape_helper(
            newshape, false, function_name_for_messages);
    }

    /**
     * This is for SOMADataFrame.  While resize_soma_joinid_shape allows the
     * user to do up the soma_joinid domain slot, without needing to specify
     * the rest (which is the common operation for experiment-level resize)
     * this allows the full-generality resize-every-index-column case
     * (which only applies to variant-indexed/non-standard dataframes).
     */
    void change_domain(
        const ArrowTable& newdomain, std::string function_name_for_messages) {
        _set_domain_helper(newdomain, true, function_name_for_messages);
    }

    /**
     * This is for SOMADataFrame.  While upgrade_soma_joinid_shape allows the
     * user to do up the soma_joinid domain slot, without needing to specify
     * the rest (which is the common operation for experiment-level resize)
     * this allows the full-generality resize-every-index-column case
     * (which only applies to variant-indexed/non-standard dataframes).
     */
    void upgrade_domain(
        const ArrowTable& newdomain, std::string function_name_for_messages) {
        _set_domain_helper(newdomain, false, function_name_for_messages);
    }

    std::shared_ptr<SOMAColumn> get_column(std::string_view name) const;

    std::shared_ptr<SOMAColumn> get_column(std::size_t index) const;

   protected:
    // See top-of-file notes regarding methods for SOMADataFrame being
    // defined in this file.
    //
    // These return the shape and maxshape slots for the soma_joinid dim, if
    // the array has one. These are important test-points and dev-internal
    // access-points, in particular, for the tiledbsoma-io experiment-level
    // resizer.
    std::optional<int64_t> _maybe_soma_joinid_shape();
    std::optional<int64_t> _maybe_soma_joinid_maxshape();

    static Array _create(
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
    std::shared_ptr<Array> arr_;

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    /**
     * The caller must check the return value for .is_empty() to see if this is
     * a new-style array with current-domain support (.is_empty() is false) , or
     * an old-style array without current-domain support (.is_empty() is true).
     * We could implement this as a std::optional<CurrentDomain> return value
     * here, but, that would be a redundant indicator.
     */
    CurrentDomain _get_current_domain() const {
        return tiledb::ArraySchemaExperimental::current_domain(
            *ctx_->tiledb_ctx(), *schema_);
    }

    /**
     * Returns the core current domain in its entirety, as an Arrow
     * table for return to Python/R.
     */
    ArrowTable _get_core_current_domain() {
        return _get_core_domainish(Domainish::kind_core_current_domain);
    }

    /**
     * Returns the core domain in its entirety, as an Arrow
     * table for return to Python/R.
     */
    ArrowTable _get_core_domain() {
        return _get_core_domainish(Domainish::kind_core_domain);
    }

    /**
     * This is a code-dedupe helper for can_resize and can_upgrade_shape.
     */
    StatusAndReason _can_set_shape_helper(
        const std::vector<int64_t>& newshape,
        bool must_already_have,
        std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for can_change_domain and
     * can_upgrade_domain.
     */
    StatusAndReason _can_set_domain_helper(
        const ArrowTable& newdomain,
        bool must_already_have,
        std::string function_name_for_messages);

    /**
     * This is a second-level code-dedupe helper for _can_set_shape_helper.
     */
    StatusAndReason _can_set_shape_domainish_subhelper(
        const std::vector<int64_t>& newshape,
        std::string function_name_for_messages);

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
        int64_t newshape,
        bool must_already_have,
        std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for resize and upgrade_shape.
     */
    void _set_shape_helper(
        const std::vector<int64_t>& newshape,
        bool must_already_have,
        std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for resize_soma_joinid_shape and
     * upgrade_soma_joinid_shape.
     */
    void _set_soma_joinid_shape_helper(
        int64_t newshape,
        bool must_already_have,
        std::string function_name_for_messages);

    /**
     * This is a code-dedupe helper method for change_domain and upgrade_domain.
     */
    void _set_domain_helper(
        const ArrowTable& newdomain,
        bool must_already_have,
        std::string function_name_for_messages);

    /**
     * This is a helper for can_upgrade_domain.
     */
    template <typename T>
    StatusAndReason _can_set_dataframe_domainish_slot_checker_non_string(
        bool check_current_domain,
        const ArrowTable& domain_table,
        std::string dim_name) {
        std::pair<T, T> old_lo_hi = check_current_domain ?
                                        _core_current_domain_slot<T>(dim_name) :
                                        _core_domain_slot<T>(dim_name);
        std::vector<T>
            new_lo_hi = ArrowAdapter::get_table_non_string_column_by_name<T>(
                domain_table, dim_name);
        if (new_lo_hi.size() != 2) {
            throw TileDBSOMAError(
                "internal coding error detected at "
                "_can_set_dataframe_domainish_slot_checker");
        }

        const T& old_lo = old_lo_hi.first;
        const T& old_hi = old_lo_hi.second;
        const T& new_lo = new_lo_hi[0];
        const T& new_hi = new_lo_hi[1];

        // If we're checking against the core current domain: the user-provided
        // domain must contain the core current domain.
        //
        // If we're checking against the core (max) domain: the user-provided
        // domain must be contained within the core (max) domain.

        if (new_lo > new_hi) {
            return std::pair(
                false,
                "index-column name " + dim_name + ": new lower > new upper");
        }

        if (check_current_domain) {
            if (new_lo > old_lo) {
                return std::pair(
                    false,
                    "index-column name " + dim_name +
                        ": new lower > old lower (downsize is unsupported)");
            }
            if (new_hi < old_hi) {
                return std::pair(
                    false,
                    "index-column name " + dim_name +
                        ": new upper < old upper (downsize is unsupported)");
            }
        } else {
            if (new_lo < old_lo) {
                return std::pair(
                    false,
                    "index-column name " + dim_name +
                        ": new lower < limit lower");
            }
            if (new_hi > old_hi) {
                return std::pair(
                    false,
                    "index-column name " + dim_name +
                        ": new upper > limit upper");
            }
        }
        return std::pair(true, "");
    }

    /**
     * This is a helper for can_upgrade_domain.
     */
    StatusAndReason _can_set_dataframe_domainish_slot_checker_string(
        bool /*check_current_domain*/,
        const ArrowTable& domain_table,
        std::string dim_name) {
        std::vector<std::string>
            new_lo_hi = ArrowAdapter::get_table_string_column_by_name(
                domain_table, dim_name);
        if (new_lo_hi.size() != 2) {
            throw TileDBSOMAError(
                "internal coding error detected at "
                "_can_set_dataframe_domainish_slot_checker");
        }

        const std::string& new_lo = new_lo_hi[0];
        const std::string& new_hi = new_lo_hi[1];

        if (new_lo != "" || new_hi != "") {
            return std::pair(
                false,
                "domain cannot be set for string index columns: please use "
                "(\"\", \"\")");
        }

        return std::pair(true, "");
    }

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

    void fill_metadata_cache(std::optional<TimestampRange> timestamp);

    void fill_columns();

    // Metadata cache
    std::map<std::string, MetadataValue> metadata_;

    // SOMAColumn list
    std::vector<std::shared_ptr<SOMAColumn>> columns_;

    // Read timestamp range (start, end)
    std::optional<TimestampRange> timestamp_;

    // The TileDB ArraySchema. The schema is inaccessible when the TileDB Array
    // is closed or opened in write mode which means we cannot use arr->schema()
    // directly in those cases. Here, we store a copy of the schema so that it
    // can be accessed in any mode
    std::shared_ptr<ArraySchema> schema_;

    // Array associated with metadata_. Metadata values need to be
    // accessible in write mode as well. We need to keep this read-mode
    // array alive in order for the metadata value pointers in the cache to
    // be accessible
    std::shared_ptr<Array> meta_cache_arr_;

    // Unoptimized method for computing nnz() (issue `count_cells` query)
    uint64_t _nnz_slow();
};

}  // namespace tiledbsoma

#endif  // SOMA_ARRAY
