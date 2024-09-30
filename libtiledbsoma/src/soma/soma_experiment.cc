/**
 * @file   soma_experiment.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
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
 *   This file defines the SOMAExperiment class.
 */

#include "soma_experiment.h"
#include "soma_collection.h"
#include "soma_measurement.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAExperiment::create(
    std::string_view uri,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    try {
        std::filesystem::path experiment_uri(uri);

        SOMAGroup::create(
            ctx, experiment_uri.string(), "SOMAExperiment", timestamp);
        SOMADataFrame::create(
            (experiment_uri / "obs").string(),
            std::move(schema),
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            platform_config,
            timestamp);
        SOMACollection::create(
            (experiment_uri / "ms").string(), ctx, timestamp);

        auto name = std::string(std::filesystem::path(uri).filename());
        auto group = SOMAGroup::open(
            OpenMode::write, experiment_uri.string(), ctx, name, timestamp);
        group->set(
            (experiment_uri / "obs").string(),
            URIType::absolute,
            "obs",
            "SOMADataFrame");
        group->set(
            (experiment_uri / "ms").string(),
            URIType::absolute,
            "ms",
            "SOMACollection");
        group->close();
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAExperiment> SOMAExperiment::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    try {
        return std::make_unique<SOMAExperiment>(mode, uri, ctx, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::shared_ptr<SOMADataFrame> SOMAExperiment::obs(
    std::vector<std::string> column_names, ResultOrder result_order) {
    if (obs_ == nullptr) {
        obs_ = SOMADataFrame::open(
            (std::filesystem::path(uri()) / "obs").string(),
            OpenMode::read,
            ctx(),
            column_names,
            result_order,
            timestamp());
    }
    return obs_;
}

std::shared_ptr<SOMACollection> SOMAExperiment::ms() {
    if (ms_ == nullptr) {
        ms_ = SOMACollection::open(
            (std::filesystem::path(uri()) / "ms").string(),
            OpenMode::read,
            ctx(),
            timestamp());
    }
    return ms_;
}

}  // namespace tiledbsoma
