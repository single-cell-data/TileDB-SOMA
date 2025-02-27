/**
 * @file   soma_experiment.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
    const std::unique_ptr<ArrowSchema>& schema,
    const ArrowTable& index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    try {
        std::filesystem::path experiment_uri(uri);

        SOMAGroup::create(
            ctx, experiment_uri.string(), "SOMAExperiment", timestamp);
        SOMADataFrame::create(
            (experiment_uri / "obs").string(),
            schema,
            index_columns,
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
        auto group = std::make_unique<SOMAExperiment>(
            mode, uri, ctx, timestamp);

        if (!group->check_type("SOMAExperiment")) {
            throw TileDBSOMAError(
                "[SOMAExperiment::open] Object is not a SOMAExperiment");
        }

        return group;
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::shared_ptr<SOMADataFrame> SOMAExperiment::obs() {
    if (obs_ == nullptr) {
        obs_ = SOMADataFrame::open(
            (std::filesystem::path(uri()) / "obs").string(),
            OpenMode::read,
            ctx(),
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
