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
//= public non-static
//===================================================================

SOMAGroup::SOMAGroup(
    std::shared_ptr<Context> ctx,
    const std::string& uri,
    tiledb_query_type_t mode)
    : ctx_(ctx) {
    group_ = std::make_unique<Group>(*ctx_, uri, mode);
}

void SOMAGroup::close() {
    group_.get()->close();
}

std::string SOMAGroup::uri() const {
    return group_.get()->uri();
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
    const std::string& uri, const bool& relative, const std::string& name) {
    group_.get()->add_member(uri, relative, name);
}

uint64_t SOMAGroup::get_length() const {
    return group_.get()->member_count();
}

void SOMAGroup::remove_member(const std::string& name) {
    group_.get()->remove_member(name);
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