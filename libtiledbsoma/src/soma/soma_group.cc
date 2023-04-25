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

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAGroup::create(std::shared_ptr<Context> ctx, const std::string& uri) {
    Group::create(*ctx, uri);
}

std::unique_ptr<SOMAGroup> SOMAGroup::open(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::map<std::string, std::string> platform_config) {
    return std::make_unique<SOMAGroup>(
        mode, uri, std::make_shared<Context>(Config(platform_config)));
}

std::unique_ptr<SOMAGroup> SOMAGroup::open(
    tiledb_query_type_t mode,
    std::shared_ptr<Context> ctx,
    std::string_view uri) {
    return std::make_unique<SOMAGroup>(mode, uri, ctx);
}

//===================================================================
//= public non-static
//===================================================================

SOMAGroup::SOMAGroup(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::shared_ptr<Context> ctx)
    : ctx_(ctx) {
    group_ = std::make_unique<Group>(*ctx_, std::string(uri), mode);
}

void SOMAGroup::open(tiledb_query_type_t query_type) {
    group_.get()->open(query_type);
}

void SOMAGroup::close() {
    group_.get()->close();
}

std::string SOMAGroup::uri() const {
    return group_.get()->uri();
}

std::shared_ptr<Context> SOMAGroup::ctx() {
    return ctx_;
}

tiledb::Object SOMAGroup::get_member(uint64_t index) const {
    return group_.get()->member(index);
}

tiledb::Object SOMAGroup::get_member(const std::string& name) const {
    return group_.get()->member(name);
}

bool SOMAGroup::has_member(const std::string& name) {
    try {
        group_.get()->member(name);
    } catch (const TileDBError& e) {
        return false;
    }
    return true;
}

void SOMAGroup::add_member(
    const std::string& uri, bool relative, const std::string& name) {
    group_.get()->add_member(uri, relative, name);
}

uint64_t SOMAGroup::get_length() const {
    return group_.get()->member_count();
}

void SOMAGroup::remove_member(const std::string& name) {
    group_.get()->remove_member(name);
}

void SOMAGroup::set_metadata(
    const std::string& key,
    tiledb_datatype_t value_type,
    uint32_t value_num,
    const void* value) {
    group_.get()->put_metadata(key, value_type, value_num, value);
}

void SOMAGroup::delete_metadata(const std::string& key) {
    group_.get()->delete_metadata(key);
}

void SOMAGroup::get_metadata(
    const std::string& key,
    tiledb_datatype_t* value_type,
    uint32_t* value_num,
    const void** value) {
    group_.get()->get_metadata(key, value_type, value_num, value);
}

void SOMAGroup::get_metadata_from_index(
    uint64_t index,
    std::string* key,
    tiledb_datatype_t* value_type,
    uint32_t* value_num,
    const void** value) {
    group_.get()->get_metadata_from_index(
        index, key, value_type, value_num, value);
}

bool SOMAGroup::has_metadata(
    const std::string& key, tiledb_datatype_t* value_type) {
    return group_.get()->has_metadata(key, value_type);
}

uint64_t SOMAGroup::metadata_num() const {
    return group_.get()->metadata_num();
}

std::map<std::string, std::string> SOMAGroup::member_to_uri_mapping() const {
    std::map<std::string, std::string> result;
    for (uint64_t i = 0; i < this->get_length(); ++i) {
        auto mem = this->get_member(i);
        result[mem.name().value()] = mem.uri();
    }
    return result;
}

}  // namespace tiledbsoma