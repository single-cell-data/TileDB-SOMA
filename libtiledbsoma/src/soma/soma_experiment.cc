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
#include "soma_dataframe.h"
#include "soma_measurement.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAExperiment::create(
    std::string_view uri, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    try {
        // Root SOMA objects include a `dataset_type` entry to allow the TileDB Cloud UI to detect that they are SOMA datasets.
        SOMAGroup::create(ctx, uri, "SOMAExperiment", {{"dataset_type", "soma"}}, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAExperiment> SOMAExperiment::open(
    std::string_view uri, OpenMode mode, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    try {
        return std::make_unique<SOMAExperiment>(mode, uri, ctx, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::shared_ptr<SOMADataFrame> SOMAExperiment::obs() {
    if (obs_ == nullptr) {
        obs_ = SOMADataFrame::open(
            (std::filesystem::path(uri()) / "obs").string(), OpenMode::soma_read, ctx(), timestamp());
    }
    return obs_;
}

std::shared_ptr<SOMACollection> SOMAExperiment::ms() {
    if (ms_ == nullptr) {
        ms_ = SOMACollection::open(
            (std::filesystem::path(uri()) / "ms").string(), OpenMode::soma_read, ctx(), timestamp());
    }
    return ms_;
}

}  // namespace tiledbsoma
