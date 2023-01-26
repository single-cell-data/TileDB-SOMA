/**
 * @file   soma_reader.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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
 *   This file defines the SOMAReader class.
 */

#include "tiledbsoma/soma_reader.h"
#include "tiledbsoma/logger_public.h"
#include "tiledbsoma/util.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMAReader> SOMAReader::open(
    std::string_view uri,
    std::string_view name,
    std::map<std::string, std::string> platform_config,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMAReader>(
        uri,
        name,
        std::make_shared<Context>(Config(platform_config)),
        column_names,
        batch_size,
        result_order,
        timestamp);
}

std::unique_ptr<SOMAReader> SOMAReader::open(
    std::shared_ptr<Context> ctx,
    std::string_view uri,
    std::string_view name,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMAReader>(
        uri, name, ctx, column_names, batch_size, result_order, timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMAReader::SOMAReader(
    std::string_view uri,
    std::string_view name,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(uri))
    , timestamp_(timestamp) {
    // Validate parameters
    try {
        LOG_DEBUG(fmt::format("[SOMAReader] opening array '{}'", uri_));
        auto array = std::make_shared<Array>(*ctx_, uri_, TILEDB_READ);
        if (timestamp) {
            if (timestamp->first > timestamp->second) {
                throw std::invalid_argument("timestamp start > end");
            }
            array->set_open_timestamp_start(timestamp->first);
            array->set_open_timestamp_end(timestamp->second);
            array->reopen();
        }
        mq_ = std::make_unique<ManagedQuery>(array, name);
        LOG_DEBUG(
            fmt::format("timestamp_start = {}", array->open_timestamp_start()));
        LOG_DEBUG(
            fmt::format("timestamp_end = {}", array->open_timestamp_end()));
    } catch (const std::exception& e) {
        throw TileDBSOMAError(
            fmt::format("Error opening array: '{}'\n  {}", uri_, e.what()));
    }

    reset(column_names, batch_size, result_order);
}

void SOMAReader::reset(
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order) {
    // Reset managed query
    mq_->reset();

    if (!column_names.empty()) {
        mq_->select_columns(column_names);
    }

    batch_size_ = batch_size;

    result_order_ = "auto";
    if (result_order != "auto") {  // default "auto" is set in soma_reader.h
        tiledb_layout_t layout;
        if (result_order == "row-major") {
            layout = TILEDB_ROW_MAJOR;
        } else if (result_order == "column-major") {
            layout = TILEDB_COL_MAJOR;
        } else {
            throw TileDBSOMAError(
                fmt::format("Unknown result_order '{}'", result_order));
        }
        mq_->set_layout(layout);
        result_order_ = result_order;
    }

    first_read_next_ = true;
    submitted_ = false;
}

void SOMAReader::submit() {
    // Submit the query
    mq_->submit();
    submitted_ = true;
}

std::optional<std::shared_ptr<ArrayBuffers>> SOMAReader::read_next() {
    if (!submitted_) {
        throw TileDBSOMAError(
            "[SOMAReader] submit must be called before read_next");
    }

    // Always return results from the first call to read_next()
    if (first_read_next_) {
        first_read_next_ = false;
        return mq_->results();
    }

    // If the query is complete, return `std::nullopt`.
    if (mq_->is_complete()) {
        return std::nullopt;
    }

    // Submit the query
    mq_->submit();

    // Return the results, possibly incomplete
    return mq_->results();
}

uint64_t SOMAReader::nnz() {
    // Verify array is sparse
    if (mq_->schema()->array_type() != TILEDB_SPARSE) {
        throw TileDBSOMAError(
            "[SOMAReader] nnz is only supported for sparse arrays");
    }

    // Load fragment info
    FragmentInfo fragment_info(*ctx_, uri_);
    fragment_info.load();

    LOG_DEBUG(fmt::format("[SOMAReader] Fragment info for array '{}'", uri_));
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
        if (frag_ts.first != frag_ts.second) {
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
            "[SOMAReader] fragment {} non-empty domain = [{}, {}]",
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
            "[SOMAReader] Checking {} < {}",
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

uint64_t SOMAReader::nnz_slow() {
    // If duplicates are allowed, we cannot count simply count cells.
    if (mq_->schema()->allows_dups()) {
        throw TileDBSOMAError(
            "[SOMAReader] nnz not supported when duplicates are allowed");
    }

    LOG_WARN(
        "[SOMAReader] nnz() found consolidated or overlapping fragments, "
        "counting cells...");

    auto sr = SOMAReader::open(
        ctx_,
        uri_,
        "count_cells",
        {mq_->schema()->domain().dimension(0).name()},
        batch_size_,
        result_order_,
        timestamp_);
    sr->submit();

    uint64_t total_cell_num = 0;
    while (auto batch = sr->read_next()) {
        total_cell_num += (*batch)->num_rows();
    }

    return total_cell_num;
}

}  // namespace tiledbsoma
