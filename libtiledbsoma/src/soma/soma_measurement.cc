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
#include <tiledb/tiledb>
#include "soma_collection.h"
#include "soma_dataframe.h"
#include "soma_experiment.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAMeasurement::create(
    std::string_view uri, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    try {
        SOMAGroup::create(ctx, uri, "SOMAMeasurement", {}, timestamp);
    } catch (tiledb::TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAMeasurement> SOMAMeasurement::open(
    std::string_view uri, OpenMode mode, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    try {
        return std::make_unique<SOMAMeasurement>(mode, uri, ctx, timestamp);
    } catch (tiledb::TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::shared_ptr<SOMADataFrame> SOMAMeasurement::var() {
    return std::dynamic_pointer_cast<SOMADataFrame>(get("var"));
    ;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::X() {
    return std::dynamic_pointer_cast<SOMACollection>(get("X"));
}

std::shared_ptr<SOMACollection> SOMAMeasurement::obsm() {
    return std::dynamic_pointer_cast<SOMACollection>(get("obsm"));
}

std::shared_ptr<SOMACollection> SOMAMeasurement::obsp() {
    return std::dynamic_pointer_cast<SOMACollection>(get("obsp"));
}

std::shared_ptr<SOMACollection> SOMAMeasurement::varm() {
    return std::dynamic_pointer_cast<SOMACollection>(get("varm"));
}

std::shared_ptr<SOMACollection> SOMAMeasurement::varp() {
    return std::dynamic_pointer_cast<SOMACollection>(get("varp"));
}
}  // namespace tiledbsoma
