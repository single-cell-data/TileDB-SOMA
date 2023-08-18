/**
 * @file   soma_group.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAGroup class.
 */

#include "soma_group.h"
#include "../soma/logger_public.h"
#include "../utils/util.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAGroup::create(
    std::shared_ptr<Context> ctx, std::string_view uri, std::string soma_type) {
    Group::create(*ctx, std::string(uri));
    auto group = Group(*ctx, std::string(uri), TILEDB_WRITE);
    group.put_metadata(
        "soma_object_type", TILEDB_STRING_UTF8, 1, soma_type.c_str());
    group.close();
}

std::unique_ptr<SOMAGroup> SOMAGroup::open(
    OpenMode mode,
    std::string_view uri,
    std::string_view name,
    std::map<std::string, std::string> platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMAGroup>(
        mode,
        uri,
        name,
        std::make_shared<Context>(Config(platform_config)),
        timestamp);
}

std::unique_ptr<SOMAGroup> SOMAGroup::open(
    OpenMode mode,
    std::shared_ptr<Context> ctx,
    std::string_view uri,
    std::string_view name,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMAGroup>(mode, uri, name, ctx, timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMAGroup::SOMAGroup(
    OpenMode mode,
    std::string_view uri,
    std::string_view name,
    std::shared_ptr<Context> ctx,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(uri))
    , name_(name) {
    auto cfg = ctx_->config();
    if (timestamp) {
        if (timestamp->first > timestamp->second) {
            throw std::invalid_argument("timestamp start > end");
        }
        cfg["sm.group.timestamp_start"] = timestamp->first;
        cfg["sm.group.timestamp_end"] = timestamp->second;
    }
    group_ = std::make_unique<Group>(
        *ctx_,
        std::string(uri),
        mode == OpenMode::read ? TILEDB_READ : TILEDB_WRITE,
        cfg);

    fill_caches();
}

void SOMAGroup::fill_caches() {
    std::shared_ptr<Group> grp;
    if (group_->query_type() == TILEDB_WRITE) {
        grp = std::make_shared<Group>(*ctx_, uri_, TILEDB_READ);
    } else {
        grp = group_;
    }

    for (uint64_t idx = 0; idx < grp->metadata_num(); ++idx) {
        std::string key;
        tiledb_datatype_t value_type;
        uint32_t value_num;
        const void* value;
        grp->get_metadata_from_index(
            idx, &key, &value_type, &value_num, &value);
        MetadataValue mdval(value_type, value_num, value);
        std::pair<std::string, const MetadataValue> mdpair(key, mdval);
        metadata_.insert(mdpair);
    }

    for (uint64_t i = 0; i < grp->member_count(); ++i) {
        auto mem = grp->member(i);
        member_to_uri_[mem.name().value()] = mem.uri();
    }

    if (group_->query_type() == TILEDB_WRITE) {
        grp->close();
    }
}

void SOMAGroup::open(
    OpenMode query_type,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto cfg = ctx_->config();
    if (timestamp) {
        if (timestamp->first > timestamp->second) {
            throw std::invalid_argument("timestamp start > end");
        }
        cfg["sm.group.timestamp_start"] = timestamp->first;
        cfg["sm.group.timestamp_end"] = timestamp->second;
    }
    group_->set_config(cfg);
    group_->open(query_type == OpenMode::read ? TILEDB_READ : TILEDB_WRITE);
}

void SOMAGroup::close() {
    group_->close();
}

const std::string SOMAGroup::uri() const {
    return group_->uri();
}

std::shared_ptr<Context> SOMAGroup::ctx() {
    return ctx_;
}

tiledb::Object SOMAGroup::get_member(uint64_t index) const {
    return group_->member(index);
}

tiledb::Object SOMAGroup::get_member(const std::string& name) const {
    return group_->member(name);
}

bool SOMAGroup::has_member(const std::string& name) {
    try {
        group_->member(name);
    } catch (const TileDBError& e) {
        return false;
    }
    return true;
}

void SOMAGroup::add_member(
    const std::string& uri, URIType uri_type, const std::string& name) {
    bool relative = uri_type == URIType::relative;
    if (uri_type == URIType::automatic) {
        relative = uri.find("://") != std::string::npos;
    }
    group_->add_member(uri, relative, name);
    member_to_uri_[name] = uri;
}

uint64_t SOMAGroup::get_length() const {
    return group_->member_count();
}

void SOMAGroup::remove_member(const std::string& name) {
    group_->remove_member(name);
}

std::map<std::string, std::string> SOMAGroup::member_to_uri_mapping() const {
    return member_to_uri_;
}

void SOMAGroup::set_metadata(
    const std::string& key,
    tiledb_datatype_t value_type,
    uint32_t value_num,
    const void* value) {
    if (key.compare("soma_object_type") == 0) {
        throw TileDBSOMAError("soma_object_type cannot be modified.");
    }

    group_->put_metadata(key, value_type, value_num, value);
    MetadataValue mdval(value_type, value_num, value);
    std::pair<std::string, const MetadataValue> mdpair(key, mdval);
    metadata_.insert(mdpair);
}

void SOMAGroup::delete_metadata(const std::string& key) {
    if (key.compare("soma_object_type") == 0) {
        throw TileDBSOMAError("soma_object_type cannot be deleted.");
    }

    group_->delete_metadata(key);
    metadata_.erase(key);
}

std::map<std::string, MetadataValue> SOMAGroup::get_metadata() {
    return metadata_;
}

std::optional<MetadataValue> SOMAGroup::get_metadata(const std::string& key) {
    if (metadata_.count(key) == 0) {
        return std::nullopt;
    }
    return metadata_[key];
}

bool SOMAGroup::has_metadata(const std::string& key) {
    return metadata_.count(key) != 0;
}

uint64_t SOMAGroup::metadata_num() const {
    return metadata_.size();
}

}  // namespace tiledbsoma