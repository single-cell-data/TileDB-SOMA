/**
 * @file   soma_collection_base.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMACollectionBase class.
 */

#include "soma_collection_base.h"
#include "common/logging/impl/logger.h"
#include "soma_collection.h"
#include "soma_dataframe.h"
#include "soma_dense_ndarray.h"
#include "soma_experiment.h"
#include "soma_measurement.h"
#include "soma_object.h"
#include "soma_sparse_ndarray.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public non-static
//===================================================================

SOMACollectionBase::SOMACollectionBase(
    OpenMode mode,
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp,
    std::optional<std::string> soma_type)
    : SOMAGroup(mode, uri, ctx, std::filesystem::path(uri).filename().string(), timestamp) {
    if (soma_type.has_value()) {  // No SOMA type checking needed.
        auto type_metadata = type();
        if (!type_metadata.has_value()) {
            throw TileDBSOMAError(
                fmt::format(
                    "Unable to open a {} at '{}'. Object is missing required '{}' metadata key.",
                    soma_type.value(),
                    std::string(uri),
                    SOMA_OBJECT_TYPE_KEY));
        }
        if (type_metadata.value() != soma_type.value()) {
            throw TileDBSOMAError(
                fmt::format(
                    "Unable to open a {} at '{}'. The object at this location is a {} not a {}.",
                    soma_type.value(),
                    std::string(uri),
                    type_metadata.value(),
                    soma_type.value()));
        }
    }
    check_encoding_version();
}

void SOMACollectionBase::close() {
    for (auto mem : children_) {
        if (mem.second->is_open()) {
            mem.second->close();
        }
    }
    SOMAGroup::close();
}

std::unique_ptr<SOMAObject> SOMACollectionBase::get(const std::string& key) {
    auto tiledb_obj = SOMAGroup::get(key);
    auto soma_obj = SOMAObject::open(tiledb_obj.uri(), OpenMode::soma_read, this->ctx(), this->timestamp());
    return soma_obj;
}

std::shared_ptr<SOMACollection> SOMACollectionBase::add_new_collection(
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
    std::shared_ptr<SOMACollection> member = SOMACollection::open(uri, OpenMode::soma_read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAGroup");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMAExperiment> SOMACollectionBase::add_new_experiment(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    const common::arrow::managed_unique_ptr<ArrowSchema>& schema,
    const common::arrow::ArrowTable& index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMAExperiment::create(uri, schema, index_columns, ctx, platform_config, timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMAExperiment> member = SOMAExperiment::open(uri, OpenMode::soma_read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAGroup");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMAMeasurement> SOMACollectionBase::add_new_measurement(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    const common::arrow::managed_unique_ptr<ArrowSchema>& schema,
    const common::arrow::ArrowTable& index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMAMeasurement::create(uri, schema, index_columns, ctx, platform_config, timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMAMeasurement> member = SOMAMeasurement::open(uri, OpenMode::soma_read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAGroup");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMADataFrame> SOMACollectionBase::add_new_dataframe(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    const common::arrow::managed_unique_ptr<ArrowSchema>& schema,
    const common::arrow::ArrowTable& index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMADataFrame::create(uri, schema, index_columns, ctx, platform_config, timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMADataFrame> member = SOMADataFrame::open(uri, OpenMode::soma_read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAArray");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMADenseNDArray> SOMACollectionBase::add_new_dense_ndarray(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    std::string_view format,
    const common::arrow::ArrowTable& index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMADenseNDArray::create(uri, format, index_columns, ctx, platform_config, timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMADenseNDArray> member = SOMADenseNDArray::open(uri, OpenMode::soma_read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAArray");
    children_[std::string(key)] = member;
    return member;
}

std::shared_ptr<SOMASparseNDArray> SOMACollectionBase::add_new_sparse_ndarray(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<SOMAContext> ctx,
    std::string_view format,
    const common::arrow::ArrowTable& index_columns,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    if (!timestamp) {
        timestamp = this->timestamp();
    }

    SOMASparseNDArray::create(uri, format, index_columns, ctx, platform_config, timestamp);

    // Note that we must return a shared_ptr to the member, instead of a
    // unique_ptr because we place the SOMA object into the `children_` cache
    // in addition to returning the SOMA object to the user.
    std::shared_ptr<SOMASparseNDArray> member = SOMASparseNDArray::open(uri, OpenMode::soma_read, ctx, timestamp);
    this->set(std::string(uri), uri_type, std::string(key), "SOMAArray");
    children_[std::string(key)] = member;
    return member;
}

}  // namespace tiledbsoma
