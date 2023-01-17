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
    std::string_view result_order) {
    return std::make_unique<SOMAReader>(
        uri,
        name,
        std::make_shared<Context>(Config(platform_config)),
        column_names,
        batch_size,
        result_order);
}

std::unique_ptr<SOMAReader> SOMAReader::open(
    std::shared_ptr<Context> ctx,
    std::string_view uri,
    std::string_view name,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order) {
    return std::make_unique<SOMAReader>(
        uri, name, ctx, column_names, batch_size, result_order);
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
    std::string_view result_order)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(uri)) {
    // Validate parameters
    try {
        LOG_DEBUG(fmt::format("[SOMAReader] opening array '{}'", uri_));
        auto array = std::make_shared<Array>(*ctx_, uri_, TILEDB_READ);
        mq_ = std::make_unique<ManagedQuery>(array, name);
        LOG_DEBUG(fmt::format("timestamp = {}", array->open_timestamp_end()));
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

    //===============================================================
    // If only one fragment, return total_cell_num
    auto fragment_count = fragment_info.fragment_num();

    if (fragment_count == 0) {
        // Array schema has been created but no data have been written
        return 0;
    }

    if (fragment_count == 1) {
        return fragment_info.total_cell_num();
    }

    //===============================================================
    // Check for overlapping fragments on the first dimension
    std::vector<std::array<uint64_t, 2>> non_empty_domains(fragment_count);

    // Get all non-empty domains for first dimension
    for (uint32_t i = 0; i < fragment_count; i++) {
        // TODO[perf]: Reading fragment info is not supported on TileDB Cloud
        // yet, but reading one fragment at a time will be slow. Is there
        // another way?
        fragment_info.get_non_empty_domain(i, 0, &non_empty_domains[i]);
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

    // If multiple fragments do not overlap, return the total_cell_num
    if (!overlap) {
        return fragment_info.total_cell_num();
    }

    //===============================================================
    // Found overlapping fragments, count cells
    if (mq_->schema()->allows_dups()) {
        throw TileDBSOMAError(
            "[SOMAReader] nnz not supported for overlapping fragments with "
            "duplicates");
    }

    LOG_WARN("[SOMAReader] Found overlapping fragments, counting cells...");

    // TODO[perf]: Reuse "this" object to read, then reset the state of "this"
    // so it could be used to read again (for TileDB Cloud)
    auto sr = SOMAReader::open(
        ctx_,
        uri_,
        "count_cells",
        {mq_->schema()->domain().dimension(0).name()});
    sr->submit();

    uint64_t total_cell_num = 0;
    while (auto batch = sr->read_next()) {
        total_cell_num += (*batch)->num_rows();
    }

    return total_cell_num;
}

}  // namespace tiledbsoma
