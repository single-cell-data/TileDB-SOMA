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

std::unique_ptr<SOMAGroup> SOMAGroup::open(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::string_view name,
    std::map<std::string, std::string> platform_config,
    std::optional<uint64_t> timestamp) {
    return std::make_unique<SOMAGroup>(
        mode,
        uri,
        name,
        std::make_shared<Context>(Config(platform_config)),
        timestamp);
}

std::unique_ptr<SOMAGroup> SOMAGroup::open(
    tiledb_query_type_t mode,
    std::shared_ptr<Context> ctx,
    std::string_view uri,
    std::string_view name,
    std::optional<uint64_t> timestamp) {
    return std::make_unique<SOMAGroup>(mode, uri, name, ctx, timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMAGroup::SOMAGroup(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::string_view name,
    std::shared_ptr<Context> ctx,
    std::optional<uint64_t> timestamp)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(uri))
    , name_(name) {
    auto cfg = ctx_->config();
    if (timestamp) {
        cfg["sm.group.timestamp_end"] = timestamp.value();
    }
    group_ = std::make_unique<Group>(*ctx_, std::string(uri), mode, cfg);
}

void SOMAGroup::open(
    tiledb_query_type_t query_type, std::optional<uint64_t> timestamp) {
    if (timestamp) {
        auto cfg = ctx_->config();
        cfg["sm.group.timestamp_end"] = timestamp.value();
        group_->set_config(cfg);
    }
    group_->open(query_type);
}

void SOMAGroup::close() {
    group_->close();
}

std::string SOMAGroup::uri() const {
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
    const std::string& uri, bool relative, const std::string& name) {
    group_->add_member(uri, relative, name);
}

uint64_t SOMAGroup::get_length() const {
    return group_->member_count();
}

void SOMAGroup::remove_member(const std::string& name) {
    group_->remove_member(name);
}

std::map<std::string, std::string> SOMAGroup::member_to_uri_mapping() const {
    std::map<std::string, std::string> result;
    for (uint64_t i = 0; i < this->get_length(); ++i) {
        auto mem = this->get_member(i);
        result[mem.name().value()] = mem.uri();
    }
    return result;
}

void SOMAGroup::set_metadata(
    const std::string& key,
    tiledb_datatype_t value_type,
    uint32_t value_num,
    const void* value) {
    group_->put_metadata(key, value_type, value_num, value);
}

void SOMAGroup::delete_metadata(const std::string& key) {
    group_->delete_metadata(key);
}

MetadataValue SOMAGroup::get_metadata(const std::string& key) const {
    tiledb_datatype_t value_type;
    uint32_t value_num;
    const void* value;
    group_->get_metadata(key, &value_type, &value_num, &value);
    return MetadataValue(key, value_type, value_num, value);
}

MetadataValue SOMAGroup::get_metadata(uint64_t index) const {
    std::string key;
    tiledb_datatype_t value_type;
    uint32_t value_num;
    const void* value;
    group_->get_metadata_from_index(
        index, &key, &value_type, &value_num, &value);
    return MetadataValue(key, value_type, value_num, value);
}

bool SOMAGroup::has_metadata(const std::string& key) {
    tiledb_datatype_t value_type;
    return group_->has_metadata(key, &value_type);
}

uint64_t SOMAGroup::metadata_num() const {
    return group_->metadata_num();
}

}  // namespace tiledbsoma