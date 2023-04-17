/**
 * @file   soma_collection.cc
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
 *   This file defines the SOMACollection class.
 */

#include "../cpp_api/soma_collection.h"
#include "soma_group.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMACollection> SOMACollection::open(
    std::shared_ptr<Context> ctx, const std::string& uri) {
    return std::make_unique<SOMACollection>(ctx, uri);
}

//===================================================================
//= public non-static
//===================================================================

SOMACollection::SOMACollection(
    std::shared_ptr<Context> ctx, const std::string& uri)
    : ctx_(ctx) {
    // std::map<std::string, std::string> platform_config) {
    group_ = std::make_unique<SOMAGroup>(ctx, uri, TILEDB_READ);
}

void SOMACollection::close() {
    group_.get()->close();
}

std::string SOMACollection::uri() const {
    return group_.get()->uri();
}

// void SOMACollection::set(
//     const std::string& name, SOMAObject& object, bool use_relative_uri) {
//     group_.add_member(name, object, use_relative_uri);
// }

bool SOMACollection::has(const std::string& name) {
    return group_.get()->has_member(name);
}

uint64_t SOMACollection::count() const {
    return group_.get()->get_length();
}

void SOMACollection::del(const std::string& name) {
    group_.get()->remove_member(name);
}

std::map<std::string, std::string> SOMACollection::member_to_uri_mapping()
    const {
    return group_.get()->member_to_uri_mapping();
};

void SOMACollection::add_new_collection(
    const std::string& name, SOMACollection& collection) {
    group_.get()->add_member(collection.uri(), false, name);
}

void SOMACollection::add_new_dataframe(
    const std::string& name, SOMADataFrame& dataframe) {
    group_.get()->add_member(dataframe.uri(), false, name);
}

void SOMACollection::add_new_dense_ndarray(
    const std::string& name, const SOMADenseNDArray& array) {
    group_.get()->add_member(array.uri(), false, name);
}

void SOMACollection::add_new_sparse_ndarray(
    const std::string& name, const SOMASparseNDArray& array) {
    group_.get()->add_member(array.uri(), false, name);
}

}  // namespace tiledbsoma
