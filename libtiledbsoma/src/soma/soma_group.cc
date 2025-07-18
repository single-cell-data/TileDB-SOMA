/**
 * @file   soma_group.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAGroup class.
 */

#include "soma_group.h"
#include "../soma/logger_public.h"
#include "../utils/logger.h"  // for fmt::format
#include "../utils/util.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMAGroup> SOMAGroup::create(
    std::shared_ptr<SOMAContext> ctx,
    std::string_view uri,
    std::string_view soma_type,
    std::optional<TimestampRange> timestamp) {
    try {
        Group::create(*ctx->tiledb_ctx(), std::string(uri));

        auto group = std::make_shared<Group>(
            *ctx->tiledb_ctx(), std::string(uri), TILEDB_WRITE, _set_timestamp(ctx, timestamp));

        group->put_metadata(
            SOMA_OBJECT_TYPE_KEY, TILEDB_STRING_UTF8, static_cast<uint32_t>(soma_type.length()), soma_type.data());

        group->put_metadata(
            ENCODING_VERSION_KEY,
            TILEDB_STRING_UTF8,
            static_cast<uint32_t>(ENCODING_VERSION_VAL.length()),
            ENCODING_VERSION_VAL.c_str());

        // Root SOMA objects include a `dataset_type` entry to allow the
        // TileDB Cloud UI to detect that they are SOMA datasets.
        if (soma_type == "SOMAExperiment") {
            std::string key = "dataset_type";
            std::string dataset_type = "soma";
            group->put_metadata(
                key, TILEDB_STRING_UTF8, static_cast<uint32_t>(dataset_type.length()), dataset_type.c_str());
        }

        return std::make_unique<SOMAGroup>(ctx, group, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

std::unique_ptr<SOMAGroup> SOMAGroup::open(
    OpenMode mode,
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::string_view name,
    std::optional<TimestampRange> timestamp) {
    try {
        return std::make_unique<SOMAGroup>(mode, uri, ctx, name, timestamp);
    } catch (TileDBError& e) {
        throw TileDBSOMAError(e.what());
    }
}

//===================================================================
//= public non-static
//===================================================================

SOMAGroup::SOMAGroup(
    OpenMode mode,
    std::string_view uri,
    std::shared_ptr<SOMAContext> ctx,
    std::string_view name,
    std::optional<TimestampRange> timestamp)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(uri))
    , name_(name)
    , timestamp_(timestamp)
    , soma_mode_(mode) {
    // Note: both OpenMode.write and OpenMode.del should be opened in
    // TILEDB_WRITE mode.
    group_ = std::make_shared<Group>(
        *ctx_->tiledb_ctx(),
        std::string(uri),
        mode == OpenMode::soma_read ? TILEDB_READ : TILEDB_WRITE,
        _set_timestamp(ctx, timestamp));
    fill_caches();
}

SOMAGroup::SOMAGroup(
    std::shared_ptr<SOMAContext> ctx, std::shared_ptr<Group> group, std::optional<TimestampRange> timestamp)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(group->uri()))
    , group_(group)
    , timestamp_(timestamp) {
    switch (group_->query_type()) {
        case TILEDB_READ:
            soma_mode_ = OpenMode::soma_read;
            break;
        case TILEDB_WRITE:
            soma_mode_ = OpenMode::soma_write;
            break;
        default: {  // Only allow read/write when constructing from TileDB.
            const char* query_type_str = nullptr;
            tiledb_query_type_to_str(group_->query_type(), &query_type_str);
            throw TileDBSOMAError(
                fmt::format(
                    "Internal error: SOMAGroup constructor does not accept a "
                    "TileDB group opened in mode '{}'. The group must be opened in "
                    "either read or write mode.",
                    query_type_str));
        }
    }
    fill_caches();
}

void SOMAGroup::fill_caches() {
    if (group_->query_type() == TILEDB_WRITE) {
        cache_group_ = std::make_shared<Group>(*ctx_->tiledb_ctx(), uri_, TILEDB_READ);
    } else {
        cache_group_ = group_;
    }

    for (uint64_t idx = 0; idx < cache_group_->metadata_num(); ++idx) {
        std::string key;
        tiledb_datatype_t value_type;
        uint32_t value_num;
        const void* value;
        cache_group_->get_metadata_from_index(idx, &key, &value_type, &value_num, &value);
        MetadataValue mdval(value_type, value_num, value);
        std::pair<std::string, const MetadataValue> mdpair(key, mdval);
        metadata_.insert(mdpair);
    }

    for (uint64_t i = 0; i < cache_group_->member_count(); ++i) {
        auto mem = cache_group_->member(i);
        std::string soma_type = util::soma_type_from_tiledb_type(mem.type());
        std::string key = mem.name().has_value() ? mem.name().value() : mem.uri();
        members_map_[key] = SOMAGroupEntry(mem.uri(), soma_type);
    }
}

void SOMAGroup::open(OpenMode mode, std::optional<TimestampRange> timestamp) {
    timestamp_ = timestamp;
    soma_mode_ = mode;
    group_->set_config(_set_timestamp(ctx_, timestamp));
    // Note: both OpenMode.write and OpenMode.del should be opened in
    // TILEDB_WRITE mode.
    group_->open(mode == OpenMode::soma_read ? TILEDB_READ : TILEDB_WRITE);
    fill_caches();
}

void SOMAGroup::close() {
    if (group_->query_type() == TILEDB_WRITE)
        cache_group_->close();
    group_->close();
    metadata_.clear();
}

const std::string SOMAGroup::uri() const {
    return group_->uri();
}

std::shared_ptr<SOMAContext> SOMAGroup::ctx() {
    return ctx_;
}

tiledb::Object SOMAGroup::get(uint64_t index) const {
    return group_->member(index);
}

tiledb::Object SOMAGroup::get(const std::string& name) const {
    return group_->member(name);
}

bool SOMAGroup::has(const std::string& name) {
    try {
        group_->member(name);
    } catch (const TileDBError& e) {
        return false;
    }
    return true;
}

void SOMAGroup::set(const std::string& uri, URIType uri_type, const std::string& name, const std::string& soma_type) {
    bool relative = uri_type == URIType::relative;
    if (uri_type == URIType::automatic) {
        relative = !((uri.find("://") != std::string::npos) || (uri.find("/") == 0));
    }
    group_->add_member(uri, relative, name);
    members_map_[name] = SOMAGroupEntry(uri, soma_type);
}

uint64_t SOMAGroup::count() const {
    return group_->member_count();
}

void SOMAGroup::del(const std::string& name) {
    group_->remove_member(name);
}

std::map<std::string, SOMAGroupEntry> SOMAGroup::members_map() const {
    return members_map_;
}

std::optional<TimestampRange> SOMAGroup::timestamp() {
    return timestamp_;
}

void SOMAGroup::set_metadata(
    const std::string& key, tiledb_datatype_t value_type, uint32_t value_num, const void* value, bool force) {
    if (!force && key.compare(SOMA_OBJECT_TYPE_KEY) == 0)
        throw TileDBSOMAError(SOMA_OBJECT_TYPE_KEY + " cannot be modified.");

    if (!force && key.compare(ENCODING_VERSION_KEY) == 0)
        throw TileDBSOMAError(ENCODING_VERSION_KEY + " cannot be modified.");

    group_->put_metadata(key, value_type, value_num, value);
    MetadataValue mdval(value_type, value_num, value);
    std::pair<std::string, const MetadataValue> mdpair(key, mdval);
    metadata_.insert(mdpair);
}

void SOMAGroup::delete_metadata(const std::string& key, bool force) {
    if (!force && key.compare(SOMA_OBJECT_TYPE_KEY) == 0) {
        throw TileDBSOMAError(SOMA_OBJECT_TYPE_KEY + " cannot be deleted.");
    }

    if (!force && key.compare(ENCODING_VERSION_KEY) == 0) {
        throw TileDBSOMAError(ENCODING_VERSION_KEY + " cannot be deleted.");
    }

    group_->delete_metadata(key);
    metadata_.erase(key);
}

std::optional<MetadataValue> SOMAGroup::get_metadata(const std::string& key) {
    if (metadata_.count(key) == 0)
        return std::nullopt;

    return metadata_[key];
}

std::map<std::string, MetadataValue> SOMAGroup::get_metadata() {
    return metadata_;
}

bool SOMAGroup::has_metadata(const std::string& key) {
    return metadata_.count(key) != 0;
}

uint64_t SOMAGroup::metadata_num() const {
    return metadata_.size();
}

Config SOMAGroup::_set_timestamp(std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    auto cfg = ctx->tiledb_ctx()->config();
    if (timestamp) {
        if (timestamp->first > timestamp->second) {
            throw std::invalid_argument("timestamp start > end");
        }
        cfg["sm.group.timestamp_start"] = timestamp->first;
        cfg["sm.group.timestamp_end"] = timestamp->second;
    }
    return cfg;
}

}  // namespace tiledbsoma
