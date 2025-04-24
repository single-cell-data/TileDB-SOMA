/** @file   soma_array.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAArray class.
 */

#include "soma_array.h"
#include <tiledb/array_experimental.h>
#include <ranges>
#include "../utils/logger.h"
#include "../utils/util.h"
#include "soma_attribute.h"
#include "soma_dimension.h"
#include "soma_geometry_column.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAArray::create(
    std::shared_ptr<SOMAContext> ctx,
    std::string_view uri,
    ArraySchema schema,
    std::string_view soma_type,
    std::optional<std::string_view> soma_schema,
    std::optional<TimestampRange> timestamp) {
    _create(ctx, uri, schema, soma_type, soma_schema, timestamp);
}

Array SOMAArray::_create(
    std::shared_ptr<SOMAContext> ctx,
    std::string_view uri,
    ArraySchema schema,
    std::string_view soma_type,
    std::optional<std::string_view> soma_schema,
    std::optional<TimestampRange> timestamp) {
    // Create TileDB array.
    Array::create(std::string(uri), schema);

    // Open TileDB array at requested time.
    auto temporal_policy = timestamp.has_value() ? TemporalPolicy(
                                                       TimestampStartEnd,
                                                       timestamp->first,
                                                       timestamp->second) :
                                                   TemporalPolicy();
    Array array{
        *ctx->tiledb_ctx(), std::string(uri), TILEDB_WRITE, temporal_policy};

    // Set SOMA metadata.
    array.put_metadata(
        SOMA_OBJECT_TYPE_KEY,
        TILEDB_STRING_UTF8,
        static_cast<uint32_t>(soma_type.length()),
        soma_type.data());
    array.put_metadata(
        ENCODING_VERSION_KEY,
        TILEDB_STRING_UTF8,
        static_cast<uint32_t>(ENCODING_VERSION_VAL.length()),
        ENCODING_VERSION_VAL.c_str());
    if (soma_schema.has_value()) {
        array.put_metadata(
            TILEDB_SOMA_SCHEMA_KEY,
            TILEDB_STRING_UTF8,
            static_cast<uint32_t>(soma_schema->length()),
            soma_schema->data());
    }
    // Return internal TileDB array.
    return array;
}

std::unique_ptr<SOMAArray> SOMAArray::open(
    OpenMode mode,
    std::string_view uri,
    std::map<std::string, std::string> platform_config,
    std::optional<TimestampRange> timestamp) {
    LOG_DEBUG(
        fmt::format("[SOMAArray] static method 'cfg' opening array '{}'", uri));
    return std::make_unique<SOMAArray>(
        mode, uri, std::make_shared<SOMAContext>(platform_config), timestamp);
}

std::unique_ptr<SOMAArray> SOMAArray::open(
    OpenMode mode,
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    LOG_DEBUG(
        fmt::format("[SOMAArray] static method 'ctx' opening array '{}'", uri));
    return std::make_unique<SOMAArray>(mode, uri, ctx, timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMAArray::SOMAArray(
    OpenMode mode,
    std::string_view uri,
    std::map<std::string, std::string> platform_config,
    std::optional<TimestampRange> timestamp)
    : uri_(util::rstrip_uri(uri))
    , timestamp_(timestamp) {
    ctx_ = std::make_shared<SOMAContext>(platform_config);
    validate(mode, timestamp);
    fill_metadata_cache(timestamp);
    fill_columns();
}

SOMAArray::SOMAArray(
    OpenMode mode,
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp)
    : uri_(util::rstrip_uri(uri))
    , ctx_(ctx)
    , timestamp_(timestamp) {
    validate(mode, timestamp);
    fill_metadata_cache(timestamp);
    fill_columns();
}

SOMAArray::SOMAArray(
    std::shared_ptr<SOMAContext> ctx,
    std::shared_ptr<Array> arr,
    std::optional<TimestampRange> timestamp)
    // Ensure protected attributes initalized first in a consistent ordering
    : uri_(util::rstrip_uri(arr->uri()))
    , ctx_(ctx)
    , arr_(arr)
    // Initialize private attributes next to control the order of destruction
    , timestamp_(timestamp)
    , schema_(std::make_shared<ArraySchema>(arr->schema())) {
    fill_metadata_cache(timestamp);
    fill_columns();
}

void SOMAArray::fill_metadata_cache(std::optional<TimestampRange> timestamp) {
    if (arr_->query_type() == TILEDB_WRITE) {
        if (timestamp) {
            meta_cache_arr_ = std::make_shared<Array>(
                *ctx_->tiledb_ctx(),
                uri_,
                TILEDB_READ,
                TemporalPolicy(
                    TimestampStartEnd, timestamp->first, timestamp->second));
        } else {
            meta_cache_arr_ = std::make_shared<Array>(
                *ctx_->tiledb_ctx(), uri_, TILEDB_READ);
        }
    } else {
        meta_cache_arr_ = arr_;
    }

    metadata_.clear();

    for (uint64_t idx = 0; idx < meta_cache_arr_->metadata_num(); ++idx) {
        std::string key;
        tiledb_datatype_t value_type;
        uint32_t value_num;
        const void* value;
        meta_cache_arr_->get_metadata_from_index(
            idx, &key, &value_type, &value_num, &value);
        MetadataValue mdval(value_type, value_num, value);
        std::pair<std::string, const MetadataValue> mdpair(key, mdval);
        metadata_.insert(mdpair);
    }
}

void SOMAArray::fill_columns() {
    // Clear columns in case of reopen
    columns_.clear();

    if (type().value_or("") == "SOMAGeometryDataFrame") {
        if (!has_metadata(TILEDB_SOMA_SCHEMA_KEY)) {
            throw TileDBSOMAError(fmt::format(
                "[SOMAArray][fill_columns] Missing required metadata key '{}' "
                "from SOMAGeometryDataFrame '{}'",
                TILEDB_SOMA_SCHEMA_KEY,
                uri()));
        }
    }

    columns_ = SOMAColumn::deserialize(*ctx_->tiledb_ctx(), *arr_, metadata_);
}

const std::string SOMAArray::uri() const {
    return uri_;
};

std::shared_ptr<SOMAContext> SOMAArray::ctx() {
    return ctx_;
};

void SOMAArray::open(OpenMode mode, std::optional<TimestampRange> timestamp) {
    timestamp_ = timestamp;

    validate(mode, timestamp);
    fill_metadata_cache(timestamp_);
    fill_columns();
}

std::unique_ptr<SOMAArray> SOMAArray::reopen(
    OpenMode mode, std::optional<TimestampRange> timestamp) {
    if (arr_->query_type() == TILEDB_READ) {
        arr_->reopen();
    } else {
        arr_->close();
        arr_->open(TILEDB_WRITE);
    }
    return std::make_unique<SOMAArray>(mode, uri_, ctx_, timestamp);
}

void SOMAArray::close() {
    if (arr_->query_type() == TILEDB_WRITE) {
        meta_cache_arr_->close();
    }

    arr_->close();
    metadata_.clear();
}

uint64_t SOMAArray::ndim() const {
    return std::count_if(
        columns_.cbegin(), columns_.cend(), [](const auto& col) {
            return col->isIndexColumn();
        });
}

std::vector<std::string> SOMAArray::dimension_names() const {
    std::vector<std::string> result;
    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        result.push_back(column->name());
    }
    return result;
}

bool SOMAArray::has_dimension_name(std::string_view name) const {
    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        if (column->name() == name) {
            return true;
        }
    }
    return false;
}

std::vector<std::string> SOMAArray::attribute_names() const {
    std::vector<std::string> result;
    for (const auto& column :
         columns_ | std::views::filter([](const auto& col) {
             return !col->isIndexColumn();
         })) {
        result.push_back(column->name());
    }
    return result;
}

void SOMAArray::consolidate_and_vacuum(std::vector<std::string> modes) {
    for (auto mode : modes) {
        auto cfg = ctx_->tiledb_ctx()->config();
        cfg["sm.consolidation.mode"] = mode;
        Array::consolidate(Context(cfg), uri_);
        Array::vacuum(Context(cfg), uri_);
    }
}

void SOMAArray::set_metadata(
    const std::string& key,
    tiledb_datatype_t value_type,
    uint32_t value_num,
    const void* value,
    bool force) {
    if (!force && key.compare(SOMA_OBJECT_TYPE_KEY) == 0)
        throw TileDBSOMAError(SOMA_OBJECT_TYPE_KEY + " cannot be modified.");

    if (!force && key.compare(ENCODING_VERSION_KEY) == 0)
        throw TileDBSOMAError(ENCODING_VERSION_KEY + " cannot be modified.");

    arr_->put_metadata(key, value_type, value_num, value);

    MetadataValue mdval(value_type, value_num, value);
    std::pair<std::string, const MetadataValue> mdpair(key, mdval);
    metadata_.insert(mdpair);
}

void SOMAArray::delete_metadata(const std::string& key, bool force) {
    if (!force && key.compare(SOMA_OBJECT_TYPE_KEY) == 0) {
        throw TileDBSOMAError(SOMA_OBJECT_TYPE_KEY + " cannot be deleted.");
    }

    if (!force && key.compare(ENCODING_VERSION_KEY) == 0) {
        throw TileDBSOMAError(ENCODING_VERSION_KEY + " cannot be deleted.");
    }

    arr_->delete_metadata(key);
    metadata_.erase(key);
}

std::optional<MetadataValue> SOMAArray::get_metadata(const std::string& key) {
    if (metadata_.count(key) == 0) {
        return std::nullopt;
    }

    return metadata_[key];
}

std::map<std::string, MetadataValue> SOMAArray::get_metadata() {
    return metadata_;
}

bool SOMAArray::has_metadata(const std::string& key) {
    return metadata_.count(key) != 0;
}

uint64_t SOMAArray::metadata_num() const {
    return metadata_.size();
}

void SOMAArray::validate(
    OpenMode mode, std::optional<TimestampRange> timestamp) {
    // Validate parameters
    auto tdb_mode = mode == OpenMode::read ? TILEDB_READ : TILEDB_WRITE;

    try {
        LOG_DEBUG(fmt::format("[SOMAArray] opening array '{}'", uri_));
        if (timestamp) {
            arr_ = std::make_shared<Array>(
                *ctx_->tiledb_ctx(),
                uri_,
                tdb_mode,
                TemporalPolicy(
                    TimestampStartEnd, timestamp->first, timestamp->second));
        } else {
            arr_ = std::make_shared<Array>(*ctx_->tiledb_ctx(), uri_, tdb_mode);
        }
        LOG_TRACE(fmt::format("[SOMAArray] loading enumerations"));
        ArrayExperimental::load_all_enumerations(
            *ctx_->tiledb_ctx(), *(arr_.get()));
        schema_ = std::make_shared<ArraySchema>(arr_->schema());
    } catch (const std::exception& e) {
        throw TileDBSOMAError(
            fmt::format("Error opening array: '{}'\n  {}", uri_, e.what()));
    }
}

std::optional<TimestampRange> SOMAArray::timestamp() {
    return timestamp_;
}

Enumeration SOMAArray::get_existing_enumeration_for_column(
    std::string column_name) {
    auto tctx = ctx_->tiledb_ctx();
    auto attribute = schema_->attribute(column_name);
    auto enumeration_name = AttributeExperimental::get_enumeration_name(
        *tctx, attribute);
    if (!enumeration_name.has_value()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMAArray::get_existing_enumeration_for_column] column_name '{}' "
            "is non-enumerated",
            column_name));
    }
    return ArrayExperimental::get_enumeration(
        *tctx, *arr_, enumeration_name.value());
}

ArrowTable SOMAArray::get_enumeration_values(
    std::vector<std::string> column_names) {
    auto ncol = column_names.size();
    auto arrow_schema = ArrowAdapter::make_arrow_schema_parent(ncol);
    auto arrow_array = ArrowAdapter::make_arrow_array_parent(ncol);
    // TBD thread-pooling opportunity -- TBD if it will be worthwhile
    for (size_t i = 0; i < ncol; i++) {
        const auto column_name = column_names[i];
        auto pair = get_enumeration_values_for_column(column_name);

        arrow_array->children[i] = pair.first;
        arrow_schema->children[i] = pair.second;
    }
    return ArrowTable(std::move(arrow_array), std::move(arrow_schema));
}

std::pair<ArrowArray*, ArrowSchema*>
SOMAArray::get_enumeration_values_for_column(std::string column_name) {
    // This will throw if column name is not found
    ArrowSchema* column_schema = arrow_schema_for_column(column_name);

    // Throw if this column is not of dictionary type
    if (column_schema->dictionary == nullptr) {
        throw TileDBSOMAError(fmt::format(
            "[get_enumeration_values_for_column] column named '{}' is not of "
            "dictionary type",
            column_name));
    }

    // Release the pointers within the ArrowSchema*:
    column_schema->release(column_schema);
    // Release the ArrowSchema* itself.
    // We are needing to remember that this was allocated by malloc via
    // (as of this writing)
    // * SOMAArray::arrow_schema_for_column
    // * SOMAAttribute::arrow_schema_slot(
    // * ArrowAdapter::arrow_schema_from_tiledb_attribute
    // Follow-up audits are tracked on
    // https://github.com/single-cell-data/TileDB-SOMA/issues/3860
    free(column_schema);

    Enumeration core_enum = get_existing_enumeration_for_column(column_name);
    auto value_dtype = core_enum.type();

    ArrowSchema* output_arrow_schema = ArrowAdapter::make_arrow_schema_child(
        column_name, value_dtype);

    ArrowArray* output_arrow_array = nullptr;

    switch (value_dtype) {
        case TILEDB_STRING_UTF8:
        case TILEDB_STRING_ASCII:
        case TILEDB_CHAR:
            output_arrow_array = ArrowAdapter::make_arrow_array_child_string(
                core_enum.as_vector<std::string>());
            break;
        case TILEDB_BOOL:
            // TileDB bools are 8-bit; Arrow bools are 1-bit.
            output_arrow_array = ArrowAdapter::make_arrow_array_child_bool(
                core_enum.as_vector<uint8_t>());
            break;
        case TILEDB_INT8:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<int8_t>());
            break;
        case TILEDB_UINT8:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<uint8_t>());
            break;
        case TILEDB_INT16:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<int16_t>());
            break;
        case TILEDB_UINT16:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<uint16_t>());
            break;
        case TILEDB_INT32:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<int32_t>());
            break;
        case TILEDB_UINT32:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<uint32_t>());
            break;
        case TILEDB_INT64:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<int64_t>());
            break;
        case TILEDB_UINT64:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<uint64_t>());
            break;
        case TILEDB_FLOAT32:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<float>());
            break;
        case TILEDB_FLOAT64:
            output_arrow_array = ArrowAdapter::make_arrow_array_child(
                core_enum.as_vector<double>());
            break;

        default:
            throw TileDBSOMAError(fmt::format(
                "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                tiledb::impl::type_to_str(value_dtype)));
    }

    return std::pair(output_arrow_array, output_arrow_schema);
}

void SOMAArray::extend_enumeration_values(
    std::map<std::string, std::pair<ArrowSchema*, ArrowArray*>> values,
    bool deduplicate) {
    auto tctx = ctx_->tiledb_ctx();
    auto mq = ManagedQuery(arr_, tctx, "extend_enumeration_values");
    ArraySchemaEvolution schema_evolution(*tctx);

    // TBD thread-pooling opportunity -- TBD if it will be worthwhile
    for (const auto& [column_name, arrow_pair] : values) {
        auto [values_schema, values_array] = arrow_pair;

        if (column_name == "") {
            throw TileDBSOMAError(
                "[SOMAArray::extend_enumeration_values] column_name is the "
                "empty string");
        }

        auto attribute = schema_->attribute(column_name);
        auto enumeration_name = AttributeExperimental::get_enumeration_name(
            *tctx, attribute);
        if (!enumeration_name.has_value()) {
            throw TileDBSOMAError(fmt::format(
                "[SOMAArray::extend_enumeration_values] column_name '{}' is "
                "non-enumerated",
                column_name));
        }
        Enumeration core_enum = get_existing_enumeration_for_column(
            column_name);

        mq._extend_enumeration(
            values_schema,
            values_array,
            column_name,
            deduplicate,
            core_enum,
            schema_evolution);
    }
    schema_evolution.array_evolve(arr_->uri());
}

// Note that ArrowTable is simply our libtiledbsoma pairing of ArrowArray and
// ArrowSchema from nanoarrow.
//
// The domainish enum simply lets us re-use code which is common across
// core domain, core current domain, and core non-empty domain.
ArrowTable SOMAArray::_get_core_domainish(enum Domainish which_kind) {
    size_t array_ndim = static_cast<size_t>(ndim());

    auto arrow_schema = ArrowAdapter::make_arrow_schema_parent(array_ndim);
    auto arrow_array = ArrowAdapter::make_arrow_array_parent(array_ndim);

    size_t child_index = 0;
    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        auto kind_domain_slot = column->arrow_domain_slot(
            *ctx_, *arr_, which_kind);

        arrow_array->children[child_index] = kind_domain_slot.first;
        arrow_schema->children[child_index] = kind_domain_slot.second;

        ++child_index;
    }

    return ArrowTable(std::move(arrow_array), std::move(arrow_schema));
}

uint64_t SOMAArray::nnz(bool raise_if_slow) {
    LOG_DEBUG(fmt::format(
        "[SOMAArray] nnz with raise_if_slow='{}'",
        (raise_if_slow ? "true" : "false")));
    // Verify array is sparse
    if (schema_->array_type() != TILEDB_SPARSE) {
        throw TileDBSOMAError(
            "[SOMAArray] nnz is only supported for sparse arrays");
    }

    // Load fragment info
    FragmentInfo fragment_info(*ctx_->tiledb_ctx(), uri_);
    fragment_info.load();

    LOG_DEBUG(fmt::format("[SOMAArray] Fragment info for array '{}'", uri_));
    if (LOG_DEBUG_ENABLED()) {
        std::ostringstream fragment_info_dump;
        fragment_info_dump << fragment_info;
        fragment_info_dump.flush();

        LOG_DEBUG(fragment_info_dump.str());
    }

    // Find the subset of fragments contained within the read timestamp range or
    // are partially within [if any]
    std::vector<uint32_t> relevant_fragments;
    for (uint32_t fid = 0; fid < fragment_info.fragment_num(); fid++) {
        auto frag_ts = fragment_info.timestamp_range(fid);
        assert(frag_ts.first <= frag_ts.second);

        if (timestamp_ && (frag_ts.first > timestamp_->second ||
                           frag_ts.second < timestamp_->first)) {
            // fragment is fully outside the read timestamp range: skip it
            continue;
        }

        // fall through: fragment is fully/partially contained within the read
        // timestamp range
        relevant_fragments.push_back(fid);
    }

    auto fragment_count = relevant_fragments.size();

    if (fragment_count == 0) {
        // No data have been written [in the read timestamp range]
        return 0;
    }

    // Check for overlapping fragments on the first dimension and
    // compute total_cell_num while going through the loop
    uint64_t total_cell_num = 0;
    std::vector<std::array<uint64_t, 2>> non_empty_domains(fragment_count);

    // The loop after this only works if dim 0 is int64 soma_joinid or
    // soma_dim_0. That's the case for _almost_ all SOMADataFrame objects, but
    // not the "variant-indexed" ones: the SOMA spec only requires
    // that soma_joinid be present as a dim or an attr. It's true for all
    // SOMASparseNDArray objects.
    auto dim = tiledb_schema()->domain().dimension(0);
    auto dim_name = dim.name();
    auto type_code = dim.type();
    if ((dim_name != "soma_joinid" && dim_name != "soma_dim_0") ||
        type_code != TILEDB_INT64) {
        LOG_DEBUG(fmt::format(
            "[SOMAArray::nnz] dim 0 (type={} name={}) isn't int64 "
            "soma_joinid or int64 soma_dim_0: using _nnz_slow",
            tiledb::impl::type_to_str(type_code),
            dim_name));
        return _nnz_slow(raise_if_slow);
    }

    for (uint32_t i = 0; i < fragment_count; i++) {
        // TODO[perf]: Reading fragment info is not supported on TileDB Cloud
        // yet, but reading one fragment at a time will be slow. Is there
        // another way?
        fragment_info.get_non_empty_domain(
            relevant_fragments[i], 0, &non_empty_domains[i]);

        LOG_DEBUG(fmt::format(
            "[SOMAArray] fragment {} non-empty domain = [{}, {}]",
            i,
            non_empty_domains[i][0],
            non_empty_domains[i][1]));
    }

    // Sort non-empty domains by the start of their ranges
    // Generate sorting permutation to apply to both fragment non empty domains
    // and ids
    std::vector<size_t> permutation(relevant_fragments.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(
        permutation.begin(), permutation.end(), [&](size_t lhs, size_t rhs) {
            return non_empty_domains[lhs][0] < non_empty_domains[rhs][0];
        });

    auto permute =
        []<typename T>(
            const std::vector<T>& data,
            const std::vector<size_t>& permutation) -> std::vector<T> {
        std::vector<T> result(data.size());
        std::transform(
            permutation.cbegin(),
            permutation.cend(),
            result.begin(),
            [&](std::size_t i) { return data[i]; });
        return result;
    };

    non_empty_domains = permute(non_empty_domains, permutation);
    relevant_fragments = permute(relevant_fragments, permutation);

    // After sorting, if the end of a non-empty domain is >= the beginning of
    // the next non-empty domain, there is an overlap

    // TODO: We only need to apply the ``_nzz_slow`` on overlapping fragments
    // Tracking issue:
    // https://github.com/single-cell-data/TileDB-SOMA/issues/3908
    std::vector<std::pair<int64_t, int64_t>> overlapping_ranges;
    std::array<uint64_t, 2> current_range = non_empty_domains[0];
    bool overlap = false;

    uint64_t cell_count = fragment_info.cell_num(relevant_fragments[0]);

    auto fragment_handler = [&](uint32_t index) {
        // Current range was overlapping with previous tiles
        // and we need to add it to the overlapping ranges.
        if (overlap) {
            overlapping_ranges.push_back(
                std::make_pair(current_range[0], current_range[1]));
        } else {
            auto frag_ts = fragment_info.timestamp_range(
                relevant_fragments[index - 1]);

            if ((timestamp_ && !(frag_ts.first >= timestamp_->first &&
                                 frag_ts.second <= timestamp_->second)) ||
                (!schema_->allows_dups() && frag_ts.first != frag_ts.second)) {
                // fragment overlaps read timestamp range, but isn't fully
                // contained within: fall back to count_cells to sort that
                //                          or
                // If any relevant fragment is a consolidated fragment, fall
                // back to counting cells, because the fragment may contain
                // duplicates. If the application is allowing duplicates (in
                // which case it's the application's job to otherwise ensure
                // uniqueness), then sum-over-fragments is the right thing
                // to do.

                overlapping_ranges.push_back(
                    std::make_pair(current_range[0], current_range[1]));
            } else {
                total_cell_num += cell_count;
            }
        }
    };

    for (uint32_t i = 1; i < fragment_count; ++i) {
        if (current_range[1] < non_empty_domains[i][0]) {
            // Fragment is non-overlapping.
            fragment_handler(i);

            current_range = non_empty_domains[i];
            cell_count = fragment_info.cell_num(relevant_fragments[i]);
            overlap = false;
        } else if (current_range[1] >= non_empty_domains[i][1]) {
            // Check if range is included completely
            overlap = true;
        } else if (current_range[1] < non_empty_domains[i][1]) {
            // Check if range is partially overlapping
            current_range[1] = non_empty_domains[i][1];
            overlap = true;
        }
    }

    // Check if the last fragment is overlapped and add to overlapping_ranges
    // or add its cell count to total otherwise
    fragment_handler(relevant_fragments.size());

    // Found relevant fragments with overlap, count cells
    if (!overlapping_ranges.empty()) {
        total_cell_num += _nnz_slow(raise_if_slow, overlapping_ranges);
    }

    return total_cell_num;
}

uint64_t SOMAArray::_nnz_slow(
    bool raise_if_slow,
    const std::vector<std::pair<int64_t, int64_t>>& ranges) {
    if (raise_if_slow) {
        throw TileDBSOMAError(
            "NNZ slow path called with 'raise_if_slow==true'");
    }

    LOG_DEBUG(
        "[SOMAArray] nnz() found consolidated or overlapping fragments, "
        "counting cells...");

    std::vector<uint64_t> count(1);
    std::shared_ptr<Array> array = arr_;

    if (!arr_->is_open() || arr_->query_type() != TILEDB_READ) {
        auto temporal_policy = timestamp_.has_value() ?
                                   TemporalPolicy(
                                       TimestampStartEnd,
                                       timestamp_->first,
                                       timestamp_->second) :
                                   TemporalPolicy();
        array = std::make_shared<Array>(
            *ctx_->tiledb_ctx(), uri_, TILEDB_READ, temporal_policy);
    }

    auto query = Query(*ctx_->tiledb_ctx(), *array);
    auto subarray = Subarray(*ctx_->tiledb_ctx(), *array);

    for (const auto& range : ranges) {
        subarray.add_range(0, range.first, range.second);
    }

    // Get the default channel
    QueryChannel default_channel = QueryExperimental::get_default_channel(
        query);

    // Apply count aggregate
    default_channel.apply_aggregate("Count", CountOperation());

    query.set_layout(TILEDB_UNORDERED)
        .set_subarray(subarray)
        .set_data_buffer("Count", count);
    query.submit();

    return count[0];
}

std::vector<int64_t> SOMAArray::shape() {
    // There are two reasons for this:
    // * Transitional, non-monolithic, phased, careful development for the
    //   new-shape feature
    // * Even after the new-shape feature is fully released, there will be old
    //   arrays on disk that were created before this feature existed.
    // So this is long-term code.
    return _get_current_domain().is_empty() ?
               _shape_via_tiledb_domain() :
               _shape_via_tiledb_current_domain();
}

std::vector<int64_t> SOMAArray::maxshape() {
    return _shape_via_tiledb_domain();
}

// This is a helper for can_upgrade_shape and can_resize, which have
// much overlap.
StatusAndReason SOMAArray::_can_set_shape_helper(
    const std::vector<int64_t>& newshape,
    bool must_already_have,
    std::string function_name_for_messages) {
    // E.g. it's an error to try to upgrade_domain or resize specifying
    // a 3-D shape on a 2-D array.
    size_t arg_ndim = newshape.size();
    size_t array_ndim = static_cast<size_t>(ndim());

    if (array_ndim != arg_ndim) {
        return std::pair(
            false,
            fmt::format(
                "{}: provided shape has ndim {}, while the array has {}",
                function_name_for_messages,
                arg_ndim,
                array_ndim));
    }

    // Enforce the semantics that tiledbsoma_upgrade_shape must be called
    // only on arrays that don't have a shape set, and resize must be called
    // only on arrays that do.
    bool has_shape = has_current_domain();
    if (must_already_have) {
        // They're trying to do resize on an array that doesn't already have a
        // shape.
        if (!has_shape) {
            return std::pair(
                false,
                fmt::format(
                    "{}: array currently has no shape: please "
                    "upgrade the array.",
                    function_name_for_messages));
        }
    } else {
        // They're trying to do upgrade_shape on an array that already has a
        // shape.
        if (has_shape) {
            return std::pair(
                false,
                fmt::format(
                    "{}: array already has a shape: please use resize",
                    function_name_for_messages));
        }
    }

    // * For old-style arrays without shape: core domain (soma maxdomain) may be
    //   small (like 100) or big (like 2 billionish).
    // * For new-style arrays with shape: core current domain (soma domain) will
    //   probably be small and core domain (soma maxdomain) will be huge.
    //
    // In either case, we need to check that the user's requested shape isn't
    // outside the core domain, which is immutable. For old-style arrays,
    //
    // if the requested shape fits in the array's core domain, it's good to go
    // as a new shape.
    // For new-style arrays, we need to additionally that the the requested
    // shape (core current domain) isn't a downsize of the current one.
    auto domain_check = _can_set_shape_domainish_subhelper(
        newshape, function_name_for_messages);
    if (!domain_check.first) {
        return domain_check;
    }

    return std::pair(true, "");
}

// This is a helper for _can_set_shape_helper: it's used for comparing
// the user's requested shape against the core current domain or core (max)
// domain.
StatusAndReason SOMAArray::_can_set_shape_domainish_subhelper(
    const std::vector<int64_t>& newshape,
    std::string function_name_for_messages) {
    std::optional<NDRectangle>
        ndrect = has_current_domain() ?
                     std::make_optional<NDRectangle>(
                         tiledb::ArraySchemaExperimental::current_domain(
                             *ctx_->tiledb_ctx(), arr_->schema())
                             .ndrectangle()) :
                     std::nullopt;

    size_t idx = 0;
    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        auto status = column->can_set_current_domain_slot(
            ndrect,
            std::vector({std::make_any<std::array<int64_t, 2>>(
                std::array<int64_t, 2>({0, newshape[idx] - 1}))}));

        if (status.first == false) {
            status.second = fmt::format(
                "[{}] {}", function_name_for_messages, status.second);

            return status;
        }

        ++idx;
    }

    return std::pair(true, "");
}

StatusAndReason SOMAArray::_can_set_soma_joinid_shape_helper(
    int64_t newshape,
    bool must_already_have,
    std::string function_name_for_messages) {
    // Fail if the array doesn't already have a shape yet (they should upgrade
    // first).
    if (!must_already_have) {
        // Upgrading an array to give it a current domain
        if (has_current_domain()) {
            return std::pair(
                false,
                fmt::format(
                    "{}: dataframe already has its domain set.",
                    function_name_for_messages));
        }

    } else {
        // Resizing an array's existing current domain
        if (!has_current_domain()) {
            return std::pair(
                false,
                fmt::format(
                    "{}: dataframe currently has no domain set.",
                    function_name_for_messages));
        }
    }

    // OK if soma_joinid isn't a dim.
    if (!has_dimension_name(SOMA_JOINID)) {
        return std::pair(true, "");
    }

    // Fail if the newshape isn't within the array's core current domain.
    if (must_already_have) {
        std::pair cur_dom_lo_hi = _core_current_domain_slot<int64_t>(
            SOMA_JOINID);
        if (newshape < cur_dom_lo_hi.second) {
            return std::pair(
                false,
                fmt::format(
                    "{}: new soma_joinid shape {} < existing shape {}",
                    function_name_for_messages,
                    newshape,
                    cur_dom_lo_hi.second + 1));
        }
    }

    // Fail if the newshape isn't within the array's core (max) domain.
    std::pair dom_lo_hi = _core_domain_slot<int64_t>(SOMA_JOINID);
    if (newshape > dom_lo_hi.second) {
        return std::pair(
            false,
            fmt::format(
                "{}: new soma_joinid shape {} > maxshape {}",
                function_name_for_messages,
                newshape,
                dom_lo_hi.second + 1));
    }

    // Sucess otherwise.
    return std::pair(true, "");
}

void SOMAArray::_set_shape_helper(
    const std::vector<int64_t>& newshape,
    bool must_already_have,
    std::string function_name_for_messages) {
    if (arr_->query_type() != TILEDB_WRITE) {
        throw TileDBSOMAError(fmt::format(
            "{} array must be opened in write mode",
            function_name_for_messages));
    }

    if (!must_already_have) {
        // Upgrading an array to install a current domain
        if (!_get_current_domain().is_empty()) {
            throw TileDBSOMAError(fmt::format(
                "{}: array must not already have a shape: please upgrade it",
                function_name_for_messages));
        }
    } else {
        // Expanding an array's current domain
        if (_get_current_domain().is_empty()) {
            throw TileDBSOMAError(fmt::format(
                "{} array must already have a shape: please upgrade it",
                function_name_for_messages));
        }
    }

    // Variant-indexed dataframes must use a separate path
    _check_dims_are_int64();

    auto tctx = ctx_->tiledb_ctx();
    ArraySchemaEvolution schema_evolution(*tctx);
    CurrentDomain new_current_domain(*tctx);

    NDRectangle ndrect(*tctx, arr_->schema().domain());

    size_t array_ndim = static_cast<size_t>(ndim());
    if (newshape.size() != array_ndim) {
        throw TileDBSOMAError(fmt::format(
            "[SOMAArray::resize]: newshape has dimension count {}; array has "
            "{} ",
            newshape.size(),
            array_ndim));
    }

    size_t idx = 0;
    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        column->set_current_domain_slot(
            ndrect, std::vector<int64_t>({0, newshape[idx] - 1}));
        ++idx;
    }

    new_current_domain.set_ndrectangle(ndrect);
    schema_evolution.expand_current_domain(new_current_domain);
    schema_evolution.array_evolve(uri_);
}

void SOMAArray::_set_soma_joinid_shape_helper(
    int64_t newshape,
    bool must_already_have,
    std::string function_name_for_messages) {
    if (arr_->query_type() != TILEDB_WRITE) {
        throw TileDBSOMAError(fmt::format(
            "{}: array must be opened in write mode",
            function_name_for_messages));
    }

    if (!must_already_have) {
        // Upgrading an array to install a current domain
        if (!_get_current_domain().is_empty()) {
            throw TileDBSOMAError(fmt::format(
                "{}: array must not already have a shape",
                function_name_for_messages));
        }
    } else {
        // Expanding an array's current domain
        if (_get_current_domain().is_empty()) {
            throw TileDBSOMAError(fmt::format(
                "{} array must already have a shape",
                function_name_for_messages));
        }
    }

    auto tctx = ctx_->tiledb_ctx();
    ArraySchemaEvolution schema_evolution(*tctx);
    CurrentDomain new_current_domain(*tctx);

    if (!must_already_have) {
        // For upgrade: copy from the full/wide/max domain except for the
        // soma_joinid restriction.

        NDRectangle ndrect(*tctx, arr_->schema().domain());
        auto soma_domain = get_soma_domain();

        for (const auto& column :
             columns_ | std::views::filter([](const auto& col) {
                 return col->isIndexColumn();
             })) {
            if (column->name() == SOMA_JOINID) {
                if (column->domain_type().value() != TILEDB_INT64) {
                    throw TileDBSOMAError(fmt::format(
                        "{}: expected soma_joinid to be of type {}; got {}",
                        function_name_for_messages,
                        tiledb::impl::type_to_str(TILEDB_INT64),
                        tiledb::impl::type_to_str(
                            column->domain_type().value())));
                }

                if (column->type() !=
                    soma_column_datatype_t::SOMA_COLUMN_DIMENSION) {
                    throw TileDBSOMAError(fmt::format(
                        "{}: expected soma_joinid type to be of type "
                        "SOMA_COLUMN_DIMENSION",
                        function_name_for_messages));
                }

                column->set_current_domain_slot(
                    ndrect, std::vector<int64_t>({0, newshape - 1}));
            } else {
                column->set_current_domain_slot(
                    ndrect,
                    ArrowAdapter::get_table_any_column_by_name<2>(
                        soma_domain, column->name(), 0));
            }
        }

        new_current_domain.set_ndrectangle(ndrect);

    } else {
        // For resize: copy from the existing current domain except for the
        // new soma_joinid value.
        CurrentDomain
            old_current_domain = ArraySchemaExperimental::current_domain(
                *tctx, arr_->schema());
        NDRectangle ndrect = old_current_domain.ndrectangle();

        for (const auto& column :
             columns_ | std::views::filter([](const auto& col) {
                 return col->isIndexColumn();
             })) {
            if (column->name() == SOMA_JOINID) {
                column->set_current_domain_slot(
                    ndrect, std::vector<int64_t>({0, newshape - 1}));
                break;
            }
        }

        new_current_domain.set_ndrectangle(ndrect);
    }

    schema_evolution.expand_current_domain(new_current_domain);
    schema_evolution.array_evolve(uri_);
}

StatusAndReason SOMAArray::_can_set_domain_helper(
    const ArrowTable& newdomain,
    bool must_already_have,
    std::string function_name_for_messages) {
    // Enforce the semantics that tiledbsoma_upgrade_domain must be called
    // only on arrays that don't have a shape set, and resize must be called
    // only on arrays that do.
    if (must_already_have) {
        if (!has_current_domain()) {
            return std::pair(
                false,
                fmt::format(
                    "{}: dataframe does not have a domain: please upgrade it",
                    function_name_for_messages));
        }
    } else {
        if (has_current_domain()) {
            return std::pair(
                false,
                fmt::format(
                    "{}: dataframe already has a domain",
                    function_name_for_messages));
        }
    }

    // * For old-style dataframe without shape: core domain (soma maxdomain)
    // may
    //   be small (like 100) or big (like 2 billionish).
    // * For new-style dataframe with shape: core current domain (soma
    // domain)
    //   will probably be small and core domain (soma maxdomain) will be
    //   huge.
    //
    // In either case, we need to check that the user's requested soma
    // domain isn't outside the core domain, which is immutable. For
    // old-style dataframes, if the requested domain fits in the array's
    // core domain, it's good to go as a new soma domain.
    // For new-style dataframes, we need to additionally that the the
    // requested soma domain (core current domain) isn't a downsize of the
    // current one.
    auto current_domain_check = _can_set_dataframe_domainish_subhelper(
        newdomain, function_name_for_messages);
    if (!current_domain_check.first) {
        return current_domain_check;
    }

    return std::pair(true, "");
}

// This is a helper for can_upgrade_domain: it's used for comparing
// the user's requested soma domain against the core current domain or core
// (max) domain.
StatusAndReason SOMAArray::_can_set_dataframe_domainish_subhelper(
    const ArrowTable& newdomain, std::string function_name_for_messages) {
    if (newdomain.second->n_children != static_cast<int64_t>(ndim())) {
        return std::pair(
            false,
            fmt::format(
                "{}: requested domain has ndim={} but the dataframe has "
                "ndim={}",
                function_name_for_messages,
                newdomain.second->n_children,
                ndim()));
    }

    if (newdomain.second->n_children != newdomain.first->n_children) {
        return std::pair(
            false,
            fmt::format(
                "{}: internal coding error", function_name_for_messages));
    }

    std::optional<NDRectangle>
        ndrect = has_current_domain() ?
                     std::make_optional<NDRectangle>(
                         tiledb::ArraySchemaExperimental::current_domain(
                             *ctx_->tiledb_ctx(), arr_->schema())
                             .ndrectangle()) :
                     std::nullopt;

    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        auto status = column->can_set_current_domain_slot(
            ndrect,
            ArrowAdapter::get_table_any_column_by_name<2>(
                newdomain, column->name(), 0));

        if (status.first == false) {
            status.second = fmt::format(
                "[{}] {}", function_name_for_messages, status.second);

            return status;
        }
    }

    return std::pair(true, "");
}

void SOMAArray::_set_domain_helper(
    const ArrowTable& newdomain,
    bool must_already_have,
    std::string function_name_for_messages) {
    if (arr_->query_type() != TILEDB_WRITE) {
        throw TileDBSOMAError(fmt::format(
            "{}: array must be opened in write mode",
            function_name_for_messages));
    }

    if (must_already_have) {
        if (!has_current_domain()) {
            throw TileDBSOMAError(fmt::format(
                "{}: dataframe does not have a domain: please upgrade it",
                function_name_for_messages));
        }
    } else {
        if (has_current_domain()) {
            throw TileDBSOMAError(fmt::format(
                "{}: dataframe already has a domain",
                function_name_for_messages));
        }
    }

    if (newdomain.second->n_children != static_cast<int64_t>(ndim())) {
        throw TileDBSOMAError(fmt::format(
            "{}: requested domain has ndim={} but the dataframe has "
            "ndim={}",
            function_name_for_messages,
            newdomain.second->n_children,
            ndim()));
    }

    if (newdomain.second->n_children != newdomain.first->n_children) {
        throw TileDBSOMAError(fmt::format(
            "{}: internal coding error", function_name_for_messages));
    }

    auto tctx = ctx_->tiledb_ctx();
    NDRectangle ndrect(*tctx, arr_->schema().domain());
    CurrentDomain new_current_domain(*tctx);
    ArraySchemaEvolution schema_evolution(*tctx);

    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        column->set_current_domain_slot(
            ndrect,
            ArrowAdapter::get_table_any_column_by_name<2>(
                newdomain, column->name(), 0));
    }

    new_current_domain.set_ndrectangle(ndrect);

    schema_evolution.expand_current_domain(new_current_domain);
    schema_evolution.array_evolve(uri_);
}

std::vector<int64_t> SOMAArray::_shape_via_tiledb_current_domain() {
    // Variant-indexed dataframes must use a separate path
    _check_dims_are_int64();

    std::vector<int64_t> result;

    auto current_domain = tiledb::ArraySchemaExperimental::current_domain(
        *ctx_->tiledb_ctx(), arr_->schema());

    if (current_domain.is_empty()) {
        throw TileDBSOMAError(
            "Internal error: current domain requested for an array which "
            "does "
            "not support it");
    }

    auto t = current_domain.type();
    if (t != TILEDB_NDRECTANGLE) {
        throw TileDBSOMAError("current_domain type is not NDRECTANGLE");
    }

    NDRectangle ndrect = current_domain.ndrectangle();

    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        auto current_domain = column->core_current_domain_slot<int64_t>(ndrect);
        result.push_back(current_domain.second - current_domain.first + 1);
    }

    return result;
}

std::vector<int64_t> SOMAArray::_shape_via_tiledb_domain() {
    // Variant-indexed dataframes must use a separate path
    _check_dims_are_int64();

    std::vector<int64_t> result;
    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        auto core_domain = column->core_domain_slot<int64_t>();
        result.push_back(core_domain.second - core_domain.first + 1);
    }

    return result;
}

std::optional<int64_t> SOMAArray::_maybe_soma_joinid_shape() {
    return _get_current_domain().is_empty() ?
               _maybe_soma_joinid_shape_via_tiledb_domain() :
               _maybe_soma_joinid_shape_via_tiledb_current_domain();
}

std::optional<int64_t> SOMAArray::_maybe_soma_joinid_maxshape() {
    return _maybe_soma_joinid_shape_via_tiledb_domain();
}

std::optional<int64_t>
SOMAArray::_maybe_soma_joinid_shape_via_tiledb_current_domain() {
    if (!has_dimension_name(SOMA_JOINID)) {
        return std::nullopt;
    }

    auto column = get_column(SOMA_JOINID);
    if (column->domain_type().value() != TILEDB_INT64) {
        throw TileDBSOMAError(fmt::format(
            "expected {} dim to be {}; got {}",
            SOMA_JOINID,
            tiledb::impl::type_to_str(TILEDB_INT64),
            tiledb::impl::type_to_str(column->domain_type().value())));
    }

    auto max = column->core_current_domain_slot<int64_t>(*ctx_, *arr_).second +
               1;

    return std::optional<int64_t>(max);
}

std::optional<int64_t> SOMAArray::_maybe_soma_joinid_shape_via_tiledb_domain() {
    if (!has_dimension_name(SOMA_JOINID)) {
        return std::nullopt;
    }

    auto column = get_column(SOMA_JOINID);
    if (column->domain_type().value() != TILEDB_INT64) {
        throw TileDBSOMAError(fmt::format(
            "expected {} dim to be {}; got {}",
            SOMA_JOINID,
            tiledb::impl::type_to_str(TILEDB_INT64),
            tiledb::impl::type_to_str(column->domain_type().value())));
    }

    auto max = column->core_domain_slot<int64_t>().second + 1;

    return std::optional<int64_t>(max);
}

bool SOMAArray::_dims_are_int64() {
    for (const auto& column :
         columns_ | std::views::filter(
                        [](const auto& col) { return col->isIndexColumn(); })) {
        if (column->type() != soma_column_datatype_t::SOMA_COLUMN_DIMENSION ||
            column->domain_type().value_or(TILEDB_ANY) != TILEDB_INT64) {
            return false;
        }
    }
    return true;
}

void SOMAArray::_check_dims_are_int64() {
    if (!_dims_are_int64()) {
        throw TileDBSOMAError(
            "[SOMAArray] internal coding error: expected all dims to be "
            "int64");
    }
}

std::shared_ptr<SOMAColumn> SOMAArray::get_column(std::string_view name) const {
    auto result = std::find_if(columns_.begin(), columns_.end(), [&](auto col) {
        return col->name() == name;
    });

    if (result == columns_.end()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMAArray] internal coding error: No column named {} found",
            name));
    }

    return *result;
}

std::shared_ptr<SOMAColumn> SOMAArray::get_column(std::size_t index) const {
    if (index >= columns_.size()) {
        throw TileDBSOMAError(fmt::format(
            "[SOMAArray] internal coding error: Column index outside of range. "
            "Requested {}, but {} exist.",
            index,
            columns_.size()));
    }

    return columns_[index];
}

}  // namespace tiledbsoma
