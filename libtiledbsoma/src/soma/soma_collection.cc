/**
 * @file   soma_collection.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMACollection class.
 */

#include "soma_collection.h"
#include "soma_group.h"

namespace tiledbsoma {
using namespace tiledb;

void SOMACollection::create(
    std::string_view uri, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    try {
        SOMAGroup::create(ctx, uri, "SOMACollection", timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMACollection> SOMACollection::open(
    std::string_view uri, OpenMode mode, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    try {
        auto group = std::make_unique<SOMACollection>(mode, uri, ctx, timestamp);

        if (!group->check_type("SOMACollection")) {
            throw TileDBSOMAError("[SOMACollection::open] Object is not a SOMACollection");
        }

        return group;
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

}  // namespace tiledbsoma
