/**
 * @file   soma_array.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022-2023 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAArray class.
 */

#include "soma_array.h"
#include <tiledb/array_experimental.h>
#include "../utils/logger.h"
#include "../utils/util.h"
namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAArray::create(
    std::shared_ptr<SOMAContext> ctx,
    std::string_view uri,
    ArraySchema schema,
    std::string soma_type,
    std::optional<TimestampRange> timestamp) {
    Array::create(std::string(uri), schema);

    std::unique_ptr<Array> array;
    if (timestamp) {
        array = std::make_unique<Array>(
            *ctx->tiledb_ctx(),
            std::string(uri),
            TILEDB_WRITE,
            TemporalPolicy(
                TimestampStartEnd, timestamp->first, timestamp->second));
    } else {
        array = std::make_unique<Array>(
            *ctx->tiledb_ctx(), std::string(uri), TILEDB_WRITE);
    }

    array->put_metadata(
        SOMA_OBJECT_TYPE_KEY,
        TILEDB_STRING_UTF8,
        static_cast<uint32_t>(soma_type.length()),
        soma_type.c_str());

    array->put_metadata(
        ENCODING_VERSION_KEY,
        TILEDB_STRING_UTF8,
        static_cast<uint32_t>(ENCODING_VERSION_VAL.length()),
        ENCODING_VERSION_VAL.c_str());

    array->close();
}

std::unique_ptr<SOMAArray> SOMAArray::open(
    OpenMode mode,
    std::string_view uri,
    std::string_view name,
    std::map<std::string, std::string> platform_config,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp) {
    LOG_DEBUG(
        fmt::format("[SOMAArray] static method 'cfg' opening array '{}'", uri));
    return std::make_unique<SOMAArray>(
        mode,
        uri,
        std::make_shared<SOMAContext>(platform_config),
        name,
        column_names,
        batch_size,
        result_order,
        timestamp);
}

std::unique_ptr<SOMAArray> SOMAArray::open(
    OpenMode mode,
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::string_view name,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp) {
    LOG_DEBUG(
        fmt::format("[SOMAArray] static method 'ctx' opening array '{}'", uri));
    return std::make_unique<SOMAArray>(
        mode,
        uri,
        ctx,
        name,
        column_names,
        batch_size,
        result_order,
        timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMAArray::SOMAArray(
    OpenMode mode,
    std::string_view uri,
    std::string_view name,
    std::map<std::string, std::string> platform_config,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp)
    : uri_(util::rstrip_uri(uri))
    , result_order_(result_order)
    , timestamp_(timestamp) {
    ctx_ = std::make_shared<SOMAContext>(platform_config);
    validate(mode, name, timestamp);
    reset(column_names, batch_size, result_order);
    fill_metadata_cache();
}

SOMAArray::SOMAArray(
    OpenMode mode,
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::string_view name,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp)
    : uri_(util::rstrip_uri(uri))
    , ctx_(ctx)
    , result_order_(result_order)
    , timestamp_(timestamp) {
    validate(mode, name, timestamp);
    reset(column_names, batch_size, result_order);
    fill_metadata_cache();
}

void SOMAArray::fill_metadata_cache() {
    if (arr_->query_type() == TILEDB_WRITE) {
        meta_cache_arr_ = std::make_shared<Array>(
            *ctx_->tiledb_ctx(),
            uri_,
            TILEDB_READ,
            TemporalPolicy(
                TimestampStartEnd, timestamp()->first, timestamp()->second));
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

const std::string SOMAArray::uri() const {
    return uri_;
};

std::shared_ptr<SOMAContext> SOMAArray::ctx() {
    return ctx_;
};

void SOMAArray::open(OpenMode mode, std::optional<TimestampRange> timestamp) {
    timestamp_ = timestamp;

    validate(mode, name_, timestamp);
    reset(column_names(), batch_size_, result_order_);
    fill_metadata_cache();
}

void SOMAArray::close() {
    if (arr_->query_type() == TILEDB_WRITE)
        meta_cache_arr_->close();

    // Close the array through the managed query to ensure any pending queries
    // are completed.
    mq_->close();
    metadata_.clear();
}

void SOMAArray::reset(
    std::vector<std::string> column_names,
    std::string_view batch_size,
    ResultOrder result_order) {
    // Reset managed query
    mq_->reset();

    if (!column_names.empty()) {
        mq_->select_columns(column_names);
    }

    switch (result_order) {
        case ResultOrder::automatic:
            if (arr_->schema().array_type() == TILEDB_SPARSE)
                mq_->set_layout(TILEDB_UNORDERED);
            else
                mq_->set_layout(TILEDB_ROW_MAJOR);
            break;
        case ResultOrder::rowmajor:
            mq_->set_layout(TILEDB_ROW_MAJOR);
            break;
        case ResultOrder::colmajor:
            mq_->set_layout(TILEDB_COL_MAJOR);
            break;
        default:
            throw std::invalid_argument(fmt::format(
                "[SOMAArray] invalid ResultOrder({}) passed",
                static_cast<int>(result_order)));
    }

    batch_size_ = batch_size;
    result_order_ = result_order;
    first_read_next_ = true;
    submitted_ = false;
}

std::optional<std::shared_ptr<ArrayBuffers>> SOMAArray::read_next() {
    // If the query is complete, return `std::nullopt`
    if (mq_->is_complete(true)) {
        return std::nullopt;
    }

    // Configure query and allocate result buffers
    mq_->setup_read();

    // Continue to submit the empty query on first read to return empty results
    if (mq_->is_empty_query()) {
        if (first_read_next_) {
            first_read_next_ = false;
            return mq_->results();
        } else {
            return std::nullopt;
        }
    }

    first_read_next_ = false;

    mq_->submit_read();

    // Return the results, possibly incomplete
    return mq_->results();
}

void SOMAArray::extend_enumeration(
    std::string_view name,
    uint64_t num_elems,
    const void* data,
    uint64_t* offsets) {
    auto enmr = ArrayExperimental::get_enumeration(
        *ctx_->tiledb_ctx(), *arr_, std::string(name));

    switch (enmr.type()) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
        case TILEDB_CHAR: {
            std::vector<uint64_t> offsets_v(
                (uint32_t*)offsets, (uint32_t*)offsets + num_elems + 1);
            std::string data_v(
                (char*)data, (char*)data + offsets_v[offsets_v.size() - 1]);
            std::vector<std::string> enums_in_write;

            for (size_t offset_idx = 0; offset_idx < offsets_v.size() - 1;
                 ++offset_idx) {
                auto beg = offsets_v[offset_idx];
                auto sz = offsets_v[offset_idx + 1] - beg;
                enums_in_write.push_back(data_v.substr(beg, sz));
            }

            std::vector<std::string> extend_values;
            auto enums_existing = enmr.as_vector<std::string>();
            for (auto enum_val : enums_in_write) {
                if (std::find(
                        enums_existing.begin(),
                        enums_existing.end(),
                        enum_val) == enums_existing.end()) {
                    extend_values.push_back(enum_val);
                }
            }

            if (extend_values.size() != 0) {
                ArraySchemaEvolution se(*ctx_->tiledb_ctx());
                se.extend_enumeration(enmr.extend(extend_values));
                se.array_evolve(uri_);
            }
            break;
        }
        case TILEDB_BOOL:
        case TILEDB_INT8: {
            SOMAArray::_extend_value_helper((int8_t*)data, num_elems, enmr);
            break;
        }
        case TILEDB_UINT8: {
            SOMAArray::_extend_value_helper((uint8_t*)data, num_elems, enmr);
            break;
        }
        case TILEDB_INT16: {
            SOMAArray::_extend_value_helper((int16_t*)data, num_elems, enmr);
            break;
        }
        case TILEDB_UINT16: {
            SOMAArray::_extend_value_helper((uint16_t*)data, num_elems, enmr);
            break;
        }
        case TILEDB_INT32: {
            SOMAArray::_extend_value_helper((int32_t*)data, num_elems, enmr);
            break;
        }
        case TILEDB_UINT32: {
            SOMAArray::_extend_value_helper((uint32_t*)data, num_elems, enmr);
            break;
        }
        case TILEDB_INT64: {
            SOMAArray::_extend_value_helper((int64_t*)data, num_elems, enmr);
            break;
        }
        case TILEDB_UINT64: {
            SOMAArray::_extend_value_helper((uint64_t*)data, num_elems, enmr);
            break;
        }
        case TILEDB_FLOAT32: {
            SOMAArray::_extend_value_helper((float*)data, num_elems, enmr);
            break;
        }
        case TILEDB_FLOAT64: {
            SOMAArray::_extend_value_helper((double*)data, num_elems, enmr);
            break;
        }
        default:
            throw TileDBSOMAError(fmt::format(
                "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                tiledb::impl::type_to_str(enmr.type())));
    }
}

void SOMAArray::set_column_data(
    std::string_view name,
    uint64_t num_elems,
    const void* data,
    uint64_t* offsets,
    uint8_t* validity) {
    if (mq_->query_type() != TILEDB_WRITE) {
        throw TileDBSOMAError("[SOMAArray] array must be opened in write mode");
    }

    // Create the array_buffer_ as necessary
    if (array_buffer_ == nullptr)
        array_buffer_ = std::make_shared<ArrayBuffers>();

    // Create a ColumnBuffer object instead of passing it in as an argument to
    // `set_column_data` because ColumnBuffer::create requires a TileDB Array
    // argument which should remain a private member of SOMAArray
    auto column = ColumnBuffer::create(arr_, name);
    column->set_data(num_elems, data, offsets, validity);

    // Keep the ColumnBuffer alive by attaching it to the ArrayBuffers class
    // member. Otherwise, the data held by the ColumnBuffer will be garbage
    // collected before it is submitted to the write query
    array_buffer_->emplace(std::string(name), column);

    mq_->set_column_data(column);
};

void SOMAArray::clear_column_data() {
    array_buffer_ = nullptr;
}

void SOMAArray::write() {
    if (mq_->query_type() != TILEDB_WRITE) {
        throw TileDBSOMAError("[SOMAArray] array must be opened in write mode");
    }
    mq_->submit_write();

    mq_->reset();
    // array_buffer_ = nullptr;
}

uint64_t SOMAArray::nnz() {
    // Verify array is sparse
    if (mq_->schema()->array_type() != TILEDB_SPARSE) {
        throw TileDBSOMAError(
            "[SOMAArray] nnz is only supported for sparse arrays");
    }

    // Load fragment info
    FragmentInfo fragment_info(*ctx_->tiledb_ctx(), uri_);
    fragment_info.load();

    LOG_DEBUG(fmt::format("[SOMAArray] Fragment info for array '{}'", uri_));
    if (LOG_DEBUG_ENABLED()) {
        fragment_info.dump();
    }

    // Find the subset of fragments contained within the read timestamp range
    // [if any]
    std::vector<uint32_t> relevant_fragments;
    for (uint32_t fid = 0; fid < fragment_info.fragment_num(); fid++) {
        auto frag_ts = fragment_info.timestamp_range(fid);
        assert(frag_ts.first <= frag_ts.second);
        if (timestamp_) {
            if (frag_ts.first > timestamp_->second ||
                frag_ts.second < timestamp_->first) {
                // fragment is fully outside the read timestamp range: skip it
                continue;
            } else if (!(frag_ts.first >= timestamp_->first &&
                         frag_ts.second <= timestamp_->second)) {
                // fragment overlaps read timestamp range, but isn't fully
                // contained within: fall back to count_cells to sort that out.
                return nnz_slow();
            }
        }
        // fall through: fragment is fully contained within the read timestamp
        // range
        relevant_fragments.push_back(fid);

        // If any relevant fragment is a consolidated fragment, fall back to
        // counting cells, because the fragment may contain duplicates.
        // If the application is allowing duplicates (in which case it's the
        // application's job to otherwise ensure uniqueness), then
        // sum-over-fragments is the right thing to do.
        if (!mq_->schema()->allows_dups() && frag_ts.first != frag_ts.second) {
            return nnz_slow();
        }
    }

    auto fragment_count = relevant_fragments.size();

    if (fragment_count == 0) {
        // No data have been written [in the read timestamp range]
        return 0;
    }

    if (fragment_count == 1) {
        // Only one fragment; return its cell_num
        return fragment_info.cell_num(relevant_fragments[0]);
    }

    // Check for overlapping fragments on the first dimension and
    // compute total_cell_num while going through the loop
    uint64_t total_cell_num = 0;
    std::vector<std::array<uint64_t, 2>> non_empty_domains(fragment_count);
    for (uint32_t i = 0; i < fragment_count; i++) {
        // TODO[perf]: Reading fragment info is not supported on TileDB Cloud
        // yet, but reading one fragment at a time will be slow. Is there
        // another way?
        total_cell_num += fragment_info.cell_num(relevant_fragments[i]);
        fragment_info.get_non_empty_domain(
            relevant_fragments[i], 0, &non_empty_domains[i]);
        LOG_DEBUG(fmt::format(
            "[SOMAArray] fragment {} non-empty domain = [{}, {}]",
            i,
            non_empty_domains[i][0],
            non_empty_domains[i][1]));
    }

    // Sort non-empty domains by the start of their ranges
    std::sort(non_empty_domains.begin(), non_empty_domains.end());

    // After sorting, if the end of a non-empty domain is >= the beginning of
    // the next non-empty domain, there is an overlap
    bool overlap = false;
    for (uint32_t i = 0; i < fragment_count - 1; i++) {
        LOG_DEBUG(fmt::format(
            "[SOMAArray] Checking {} < {}",
            non_empty_domains[i][1],
            non_empty_domains[i + 1][0]));
        if (non_empty_domains[i][1] >= non_empty_domains[i + 1][0]) {
            overlap = true;
            break;
        }
    }

    // If relevant fragments do not overlap, return the total cell_num
    if (!overlap) {
        return total_cell_num;
    }
    // Found relevant fragments with overlap, count cells
    return nnz_slow();
}

uint64_t SOMAArray::nnz_slow() {
    LOG_DEBUG(
        "[SOMAArray] nnz() found consolidated or overlapping fragments, "
        "counting cells...");

    auto sr = SOMAArray::open(
        OpenMode::read,
        uri_,
        ctx_,
        "count_cells",
        {mq_->schema()->domain().dimension(0).name()},
        batch_size_,
        result_order_,
        timestamp_);

    uint64_t total_cell_num = 0;
    while (auto batch = sr->read_next()) {
        total_cell_num += (*batch)->num_rows();
    }

    return total_cell_num;
}

std::vector<int64_t> SOMAArray::shape() {
    std::vector<int64_t> result;
    auto dimensions = mq_->schema()->domain().dimensions();

    for (const auto& dim : dimensions) {
        switch (dim.type()) {
            case TILEDB_UINT8:
                result.push_back(
                    dim.domain<uint8_t>().second - dim.domain<uint8_t>().first +
                    1);
                break;
            case TILEDB_INT8:
                result.push_back(
                    dim.domain<int8_t>().second - dim.domain<int8_t>().first +
                    1);
                break;
            case TILEDB_UINT16:
                result.push_back(
                    dim.domain<uint16_t>().second -
                    dim.domain<uint16_t>().first + 1);
                break;
            case TILEDB_INT16:
                result.push_back(
                    dim.domain<int16_t>().second - dim.domain<int16_t>().first +
                    1);
                break;
            case TILEDB_UINT32:
                result.push_back(
                    dim.domain<uint32_t>().second -
                    dim.domain<uint32_t>().first + 1);
                break;
            case TILEDB_INT32:
                result.push_back(
                    dim.domain<int32_t>().second - dim.domain<int32_t>().first +
                    1);
                break;
            case TILEDB_UINT64:
                result.push_back(
                    dim.domain<uint64_t>().second -
                    dim.domain<uint64_t>().first + 1);
                break;
            case TILEDB_INT64:
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
            case TILEDB_TIME_HR:
            case TILEDB_TIME_MIN:
            case TILEDB_TIME_SEC:
            case TILEDB_TIME_MS:
            case TILEDB_TIME_US:
            case TILEDB_TIME_NS:
            case TILEDB_TIME_PS:
            case TILEDB_TIME_FS:
            case TILEDB_TIME_AS:
                result.push_back(
                    dim.domain<int64_t>().second - dim.domain<int64_t>().first +
                    1);
                break;
            default:
                throw TileDBSOMAError("Dimension must be integer type.");
        }
    }

    return result;
}

uint64_t SOMAArray::ndim() const {
    return tiledb_schema()->domain().ndim();
}

std::vector<std::string> SOMAArray::dimension_names() const {
    std::vector<std::string> result;
    auto dimensions = tiledb_schema()->domain().dimensions();
    for (const auto& dim : dimensions) {
        result.push_back(dim.name());
    }
    return result;
}

std::map<std::string, Enumeration> SOMAArray::get_attr_to_enum_mapping() {
    std::map<std::string, Enumeration> result;
    for (uint32_t i = 0; i < arr_->schema().attribute_num(); ++i) {
        auto attr = arr_->schema().attribute(i);
        if (attr_has_enum(attr.name())) {
            auto enmr_label = *get_enum_label_on_attr(attr.name());
            auto enmr = ArrayExperimental::get_enumeration(
                *ctx_->tiledb_ctx(), *arr_, enmr_label);
            result.insert({attr.name(), enmr});
        }
    }
    return result;
}

std::optional<std::string> SOMAArray::get_enum_label_on_attr(
    std::string attr_name) {
    auto attr = arr_->schema().attribute(attr_name);
    return AttributeExperimental::get_enumeration_name(
        *ctx_->tiledb_ctx(), attr);
}

bool SOMAArray::attr_has_enum(std::string attr_name) {
    return get_enum_label_on_attr(attr_name).has_value();
}

void SOMAArray::set_metadata(
    const std::string& key,
    tiledb_datatype_t value_type,
    uint32_t value_num,
    const void* value) {
    if (key.compare(SOMA_OBJECT_TYPE_KEY) == 0)
        throw TileDBSOMAError(SOMA_OBJECT_TYPE_KEY + " cannot be modified.");

    if (key.compare(ENCODING_VERSION_KEY) == 0)
        throw TileDBSOMAError(ENCODING_VERSION_KEY + " cannot be modified.");

    arr_->put_metadata(key, value_type, value_num, value);

    MetadataValue mdval(value_type, value_num, value);
    std::pair<std::string, const MetadataValue> mdpair(key, mdval);
    metadata_.insert(mdpair);
}

void SOMAArray::delete_metadata(const std::string& key) {
    if (key.compare(SOMA_OBJECT_TYPE_KEY) == 0)
        throw TileDBSOMAError(SOMA_OBJECT_TYPE_KEY + " cannot be deleted.");

    if (key.compare(ENCODING_VERSION_KEY) == 0)
        throw TileDBSOMAError(ENCODING_VERSION_KEY + " cannot be deleted.");

    arr_->delete_metadata(key);
    metadata_.erase(key);
}

std::optional<MetadataValue> SOMAArray::get_metadata(const std::string& key) {
    if (metadata_.count(key) == 0)
        return std::nullopt;

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
    OpenMode mode,
    std::string_view name,
    std::optional<TimestampRange> timestamp) {
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
        mq_ = std::make_unique<ManagedQuery>(arr_, ctx_->tiledb_ctx(), name);
    } catch (const std::exception& e) {
        throw TileDBSOMAError(
            fmt::format("Error opening array: '{}'\n  {}", uri_, e.what()));
    }
}

std::optional<TimestampRange> SOMAArray::timestamp() {
    return timestamp_;
}

}  // namespace tiledbsoma
