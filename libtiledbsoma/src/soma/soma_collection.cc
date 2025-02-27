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
#include "soma_experiment.h"
#include "soma_measurement.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMACollection::create(
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    try {
        SOMAGroup::create(ctx, uri, "SOMACollection", timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMACollection> SOMACollection::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    try {
        auto group = std::make_unique<SOMACollection>(
            mode, uri, ctx, timestamp);

        if (!group->check_type("SOMACollection")) {
            throw TileDBSOMAError(
                "[SOMACollection::open] Object is not a SOMACollection");
        }

        return group;
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

//===================================================================
//= public non-static
//===================================================================

void SOMACollection::close() {
    for (auto mem : children_) {
        if (mem.second->is_open()) {
            mem.second->close();
        }
    }
    SOMAGroup::close();
}

std::unique_ptr<SOMAObject> SOMACollection::get(const std::string& key) {
    auto tiledb_obj = SOMAGroup::get(key);
    auto soma_obj = SOMAObject::open(
        tiledb_obj.uri(), OpenMode::read, this->ctx(), this->timestamp());
    return soma_obj;
}

std::shared_ptr<SOMACollection> SOMACollection::add_new_collection(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMACollection::create(uri, ctx, timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMACollection> member = SOMACollection::open(
        uri, OpenMode::read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAGroup");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMAExperiment> SOMACollection::add_new_experiment(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMAExperiment::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        platform_config,
        timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMAExperiment> member = SOMAExperiment::open(
        uri, OpenMode::read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAGroup");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMAMeasurement> SOMACollection::add_new_measurement(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMAMeasurement::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        platform_config,
        timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMAMeasurement> member = SOMAMeasurement::open(
        uri, OpenMode::read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAGroup");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMADataFrame> SOMACollection::add_new_dataframe(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMADataFrame::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        platform_config,
        timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMADataFrame> member = SOMADataFrame::open(
        uri, OpenMode::read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAArray");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMADenseNDArray> SOMACollection::add_new_dense_ndarray(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    std::string_view format,
    ArrowTable index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMADenseNDArray::create(
        uri,
        format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        platform_config,
        timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMADenseNDArray> member = SOMADenseNDArray::open(
        uri, OpenMode::read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAArray");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMASparseNDArray> SOMACollection::add_new_sparse_ndarray(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    std::string_view format,
    ArrowTable index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMASparseNDArray::create(
        uri,
        format,
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        platform_config,
        timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMASparseNDArray> member = SOMASparseNDArray::open(
        uri, OpenMode::read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAArray");
    children_[std::string(key)] = member;
    return member;
}

}  // namespace tiledbsoma
