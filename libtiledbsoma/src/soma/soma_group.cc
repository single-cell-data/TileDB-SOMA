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
#include "../utils/util.h"
#include "common/datatype/datatype.h"
#include "common/datatype/utils.h"
#include "common/logging/impl/logger.h"
#include "common/logging/logger.h"
#include "common/metadata/metadata.h"
#include "common/metadata/utils.h"
#include "enums.h"

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

#include <filesystem>

namespace tiledbsoma {
using namespace common::type;

//==================================================================
// helper functions
//==================================================================

std::map<std::string, SOMAGroupEntry> create_member_cache(tiledb::Group& group) {
    auto get_object_type_string = [](tiledb::Object& group_member) {
        switch (group_member.type()) {
            case tiledb::Object::Type::Array:
                return "SOMAArray";
            case tiledb::Object::Type::Group:
                return "SOMAGroup";
            default:
                throw TileDBSOMAError(
                    fmt::format(
                        "Internal error: Failed to open SOMA object. Unable to resolve the TileDB object type of "
                        "member '{}'.",
                        group_member.to_str()));
        }
    };

    std::map<std::string, SOMAGroupEntry> member_cache{};
    for (uint64_t i = 0; i < group.member_count(); ++i) {
        auto mem = group.member(i);
        std::string soma_type = get_object_type_string(mem);
        std::string key = mem.name().has_value() ? mem.name().value() : mem.uri();
        member_cache[key] = SOMAGroupEntry(mem.uri(), soma_type);
    }
    return member_cache;
}

//===================================================================
//= public static
//===================================================================

void SOMAGroup::create(
    std::shared_ptr<SOMAContext> ctx,
    std::string_view uri,
    std::string_view soma_type,
    const std::unordered_map<std::string, std::string>& schema_metadata,
    std::optional<TimestampRange> timestamp) {
    ctx->validate_create_uri(uri);
    try {
        tiledb::Group::create(*ctx->tiledb_ctx(), std::string(uri));
        auto group = std::make_shared<tiledb::Group>(
            *ctx->tiledb_ctx(), std::string(uri), TILEDB_WRITE, _set_timestamp(ctx, timestamp));
        group->put_metadata(
            SOMA_OBJECT_TYPE_KEY, TILEDB_STRING_UTF8, static_cast<uint32_t>(soma_type.length()), soma_type.data());
        group->put_metadata(
            ENCODING_VERSION_KEY,
            TILEDB_STRING_UTF8,
            static_cast<uint32_t>(ENCODING_VERSION_VAL.length()),
            ENCODING_VERSION_VAL.c_str());
        for (const auto& [key, value] : schema_metadata) {
            group->put_metadata(key, TILEDB_STRING_UTF8, static_cast<uint32_t>(value.length()), value.data());
        }
        group->close();
    } catch (tiledb::TileDBError& e) {
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
    } catch (tiledb::TileDBError& e) {
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
    group_ = std::make_shared<tiledb::Group>(
        *ctx_->tiledb_ctx(),
        std::string(uri),
        mode == OpenMode::soma_read ? TILEDB_READ : TILEDB_WRITE,
        _set_timestamp(ctx, timestamp));
    cache_group_ = (group_->query_type() == TILEDB_READ) ?
                       group_ :
                       std::make_shared<tiledb::Group>(*ctx_->tiledb_ctx(), uri_, TILEDB_READ);
    metadata_cache_ = std::make_shared<common::MetadataCache>(*cache_group_);
    members_map_ = create_member_cache(*cache_group_);
}

SOMAGroup::SOMAGroup(
    std::shared_ptr<SOMAContext> ctx, std::shared_ptr<tiledb::Group> group, std::optional<TimestampRange> timestamp)
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
    cache_group_ = (group_->query_type() == TILEDB_READ) ?
                       group_ :
                       std::make_shared<tiledb::Group>(*ctx_->tiledb_ctx(), uri_, TILEDB_READ);
    metadata_cache_ = std::make_shared<common::MetadataCache>(*cache_group_);
    members_map_ = create_member_cache(*cache_group_);
}

void SOMAGroup::open(OpenMode mode, std::optional<TimestampRange> timestamp) {
    timestamp_ = timestamp;
    soma_mode_ = mode;
    group_->set_config(_set_timestamp(ctx_, timestamp));
    // Note: both OpenMode.write and OpenMode.del should be opened in
    // TILEDB_WRITE mode.
    group_->open(mode == OpenMode::soma_read ? TILEDB_READ : TILEDB_WRITE);
    cache_group_ = (group_->query_type() == TILEDB_READ) ?
                       group_ :
                       std::make_shared<tiledb::Group>(*ctx_->tiledb_ctx(), uri_, TILEDB_READ);
    metadata_cache_ = std::make_shared<common::MetadataCache>(*cache_group_);
    members_map_ = create_member_cache(*cache_group_);
    mutated_members_.clear();
}

void SOMAGroup::reopen(OpenMode mode, std::optional<TimestampRange> timestamp) {
    close(true);
    open(mode, timestamp);
}

void SOMAGroup::close([[maybe_unused]] bool recursive) {
    if (metadata_cache_)
        metadata_cache_->write(*group_);
    if (cache_group_)
        cache_group_->close();
    if (group_)
        group_->close();
}

const std::string SOMAGroup::uri() const {
    return group_->uri();
}

std::shared_ptr<SOMAContext> SOMAGroup::ctx() const {
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
    } catch (const tiledb::TileDBError& e) {
        return members_map_.contains(name);
    }
    return true;
}

void SOMAGroup::set(
    const std::string& uri,
    URIType uri_type,
    const std::string& name,
    const std::string& soma_type,
    const std::string& absolute_uri) {
    if (mutated_members_.contains(name) || members_map_.contains(name)) {
        throw std::range_error(fmt::format("replacing key '{}' is unsupported", name));
    }

    auto protocol = ctx_->tiledb_ctx()->data_protocol(uri);
    if (protocol == tiledb::Context::DataProtocol::v2) {
        auto tiledb_type = this->tiledb_type_from_soma_type(soma_type);
        bool relative = uri_type == URIType::relative;
        if (uri_type == URIType::automatic) {
            relative = !((uri.find("://") != std::string::npos) || (uri.find("/") == 0));
        }

        group_->add_member(uri, relative, name, tiledb_type == ObjectType::array ? TILEDB_ARRAY : TILEDB_GROUP);
    }

    members_map_[name] = SOMAGroupEntry(absolute_uri, soma_type);
    mutated_members_.emplace(name);
}

uint64_t SOMAGroup::count() const {
    return members_map_.size();
}

void SOMAGroup::del(const std::string& name) {
    if (mutated_members_.contains(name)) {
        throw std::range_error(fmt::format("Cannot delete previously-mutated key '{}'.", name));
    }

    group_->remove_member(name);
    members_map_.erase(name);
    mutated_members_.emplace(name);
}

std::map<std::string, SOMAGroupEntry> SOMAGroup::members_map() const {
    return members_map_;
}

std::optional<TimestampRange> SOMAGroup::timestamp() {
    return timestamp_;
}

void SOMAGroup::set_metadata(
    const std::string& key, common::DataType value_type, uint32_t value_num, const void* value, bool force) {
    if (!force && key.compare(SOMA_OBJECT_TYPE_KEY) == 0)
        throw TileDBSOMAError(SOMA_OBJECT_TYPE_KEY + " cannot be modified.");

    if (!force && key.compare(ENCODING_VERSION_KEY) == 0)
        throw TileDBSOMAError(ENCODING_VERSION_KEY + " cannot be modified.");

    metadata_cache_->set(key, common::decode_metadata(value_type, value_num, value));
}

void SOMAGroup::delete_metadata(const std::string& key, bool force) {
    if (!force && key.compare(SOMA_OBJECT_TYPE_KEY) == 0) {
        throw TileDBSOMAError(SOMA_OBJECT_TYPE_KEY + " cannot be deleted.");
    }

    if (!force && key.compare(ENCODING_VERSION_KEY) == 0) {
        throw TileDBSOMAError(ENCODING_VERSION_KEY + " cannot be deleted.");
    }

    metadata_cache_->del(key);
}

std::optional<common::MetadataValue> SOMAGroup::get_metadata(const std::string& key) {
    return metadata_cache_->get(key);
}

std::map<std::string, common::MetadataValue> SOMAGroup::get_metadata() {
    return metadata_cache_->get();
}

bool SOMAGroup::has_metadata(const std::string& key) const {
    return metadata_cache_->contains(key);
}

uint64_t SOMAGroup::metadata_num() const {
    return metadata_cache_->size();
}

std::string SOMAGroup::classname() const {
    return "Group";
}

tiledb::Config SOMAGroup::_set_timestamp(std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
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
