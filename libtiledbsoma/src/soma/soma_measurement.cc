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
    std::call_once(*var_flag_, [&]() { var_ = std::dynamic_pointer_cast<SOMADataFrame>(get_member("var")); });

    return var_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::X() {
    std::call_once(*X_flag_, [&]() { X_ = std::dynamic_pointer_cast<SOMACollection>(get_member("X")); });

    return X_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::obsm() {
    std::call_once(*obsm_flag_, [&]() { obsm_ = std::dynamic_pointer_cast<SOMACollection>(get_member("obsm")); });

    return obsm_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::obsp() {
    std::call_once(*obsp_flag_, [&]() { obsp_ = std::dynamic_pointer_cast<SOMACollection>(get_member("obsp")); });

    return obsp_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::varm() {
    std::call_once(*varm_flag_, [&]() { varm_ = std::dynamic_pointer_cast<SOMACollection>(get_member("varm")); });

    return varm_;
}

std::shared_ptr<SOMACollection> SOMAMeasurement::varp() {
    std::call_once(*varp_flag_, [&]() { varp_ = std::dynamic_pointer_cast<SOMACollection>(get_member("varp")); });

    return varp_;
}
}  // namespace tiledbsoma
