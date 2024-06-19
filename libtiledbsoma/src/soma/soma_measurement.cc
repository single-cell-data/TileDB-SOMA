/**
 * @file   soma_measurement.cc
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
 *   This file defines the SOMAMeasurement class.
 */

#include "soma_measurement.h"
#include "soma_collection.h"
#include "soma_experiment.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAMeasurement::create(
    std::string_view uri,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    try {
        std::filesystem::path measurement_uri(uri);

        SOMAGroup::create(
            ctx, measurement_uri.string(), "SOMAMeasurement", timestamp);
        SOMADataFrame::create(
            (measurement_uri / "var").string(),
            std::move(schema),
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)),
            ctx,
            platform_config,
            timestamp);
        SOMACollection::create(
            (measurement_uri / "X").string(), ctx, timestamp);
        SOMACollection::create(
            (measurement_uri / "obsm").string(), ctx, timestamp);
        SOMACollection::create(
            (measurement_uri / "obsp").string(), ctx, timestamp);
        SOMACollection::create(
            (measurement_uri / "varm").string(), ctx, timestamp);
        SOMACollection::create(
            (measurement_uri / "varp").string(), ctx, timestamp);

        auto name = std::string(std::filesystem::path(uri).filename());
        auto group = SOMAGroup::open(
            OpenMode::write, uri, ctx, name, timestamp);
        group->set(
            (measurement_uri / "var").string(),
            URIType::absolute,
            "var",
            "SOMADataFrame");
        group->set(
            (measurement_uri / "X").string(),
            URIType::absolute,
            "X",
            "SOMACollection");
        group->set(
            (measurement_uri / "obsm").string(),
            URIType::absolute,
            "obsm",
            "SOMACollection");
        group->set(
            (measurement_uri / "obsp").string(),
            URIType::absolute,
            "obsp",
            "SOMACollection");
        group->set(
            (measurement_uri / "varm").string(),
            URIType::absolute,
            "varm",
            "SOMACollection");
        group->set(
            (measurement_uri / "varp").string(),
            URIType::absolute,
            "varp",
            "SOMACollection");
        group->close();
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAMeasurement> SOMAMeasurement::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    try {
        return std::make_unique<SOMAMeasurement>(mode, uri, ctx, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::shared_ptr<SOMADataFrame> SOMAMeasurement::var(
    std::vector<std::string> column_names, ResultOrder result_order) {
    if (var_ == nullptr) {
        var_ = SOMADataFrame::open(
            (std::filesystem::path(uri()) / "var").string(),
            OpenMode::read,
            ctx(),
            column_names,
            result_order,
            timestamp());
    }
    return var_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::X() {
    if (X_ == nullptr) {
        X_ = SOMACollection::open(
            (std::filesystem::path(uri()) / "X").string(),
            OpenMode::read,
            ctx(),
            timestamp());
    }
    return X_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::obsm() {
    if (obsm_ == nullptr) {
        obsm_ = SOMACollection::open(
            (std::filesystem::path(uri()) / "obsm").string(),
            OpenMode::read,
            ctx(),
            timestamp());
    }
    return obsm_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::obsp() {
    if (obsp_ == nullptr) {
        obsp_ = SOMACollection::open(
            (std::filesystem::path(uri()) / "obsp").string(),
            OpenMode::read,
            ctx(),
            timestamp());
    }
    return obsp_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::varm() {
    if (varm_ == nullptr) {
        varm_ = SOMACollection::open(
            (std::filesystem::path(uri()) / "varm").string(),
            OpenMode::read,
            ctx(),
            timestamp());
    }
    return varm_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::varp() {
    if (varp_ == nullptr) {
        varp_ = SOMACollection::open(
            (std::filesystem::path(uri()) / "varp").string(),
            OpenMode::read,
            ctx(),
            timestamp());
    }
    return varp_;
}
}  // namespace tiledbsoma
