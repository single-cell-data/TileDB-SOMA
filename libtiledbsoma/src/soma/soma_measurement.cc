/**
 * @file   soma_measurement.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
    const std::unique_ptr<ArrowSchema>& schema,
    const ArrowTable& index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    try {
        std::filesystem::path measurement_uri(uri);

        SOMAGroup::create(
            ctx, measurement_uri.string(), "SOMAMeasurement", timestamp);
        SOMADataFrame::create(
            (measurement_uri / "var").string(),
            schema,
            index_columns,
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
        auto group = std::make_unique<SOMAMeasurement>(
            mode, uri, ctx, timestamp);

        if (!group->check_type("SOMAMeasurement")) {
            throw TileDBSOMAError(
                "[SOMAMeasurement::open] Object is not a SOMAMeasurement");
        }

        return group;
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::shared_ptr<SOMADataFrame> SOMAMeasurement::var() {
    if (var_ == nullptr) {
        var_ = SOMADataFrame::open(
            (std::filesystem::path(uri()) / "var").string(),
            OpenMode::read,
            ctx(),
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
