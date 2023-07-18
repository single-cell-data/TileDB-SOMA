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

#include "soma_collection.h"
#include "soma_experiment.h"
#include "soma_measurement.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMACollection> SOMACollection::create(
    std::shared_ptr<Context> ctx, std::string_view uri) {
    SOMAGroup::create(ctx, uri, "SOMACollection");
    return std::make_unique<SOMACollection>(TILEDB_READ, uri, ctx);
}

std::unique_ptr<SOMACollection> SOMACollection::open(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::map<std::string, std::string> platform_config) {
    return std::make_unique<SOMACollection>(
        mode, uri, std::make_shared<Context>(Config(platform_config)));
}

std::unique_ptr<SOMACollection> SOMACollection::open(
    tiledb_query_type_t mode,
    std::shared_ptr<Context> ctx,
    std::string_view uri) {
    return std::make_unique<SOMACollection>(mode, uri, ctx);
}

//===================================================================
//= public non-static
//===================================================================

SOMACollection::SOMACollection(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::shared_ptr<Context> ctx,
    std::optional<uint64_t> timestamp) {
    group_ = std::make_shared<SOMAGroup>(mode, uri, "", ctx, timestamp);
}

void SOMACollection::open(tiledb_query_type_t mode) {
    group_->open(mode);
}

void SOMACollection::close() {
    group_->close();
}

const std::string SOMACollection::uri() const {
    return group_->uri();
}

std::shared_ptr<Context> SOMACollection::ctx() {
    return group_->ctx();
}

void SOMACollection::set(
    std::string_view uri, bool relative, const std::string& key) {
    group_->add_member(std::string(uri), relative, key);
}

std::unique_ptr<SOMAObject> SOMACollection::get(const std::string& key) {
    auto member = group_->get_member(key);
    std::string soma_object_type = this->type();

    if (soma_object_type.compare("SOMACollection") == 0)
        return SOMACollection::open(TILEDB_READ, member.uri());
    else if (soma_object_type.compare("SOMAExperiment") == 0)
        return SOMAExperiment::open(TILEDB_READ, member.uri());
    else if (soma_object_type.compare("SOMAMeasurement") == 0)
        return SOMAMeasurement::open(TILEDB_READ, member.uri());
    else if (soma_object_type.compare("SOMADataFrame") == 0)
        return SOMADataFrame::open(TILEDB_READ, member.uri());
    else if (soma_object_type.compare("SOMASparseNDArray") == 0)
        return SOMASparseNDArray::open(TILEDB_READ, member.uri());
    else if (soma_object_type.compare("SOMADenseNDArray") == 0)
        return SOMADenseNDArray::open(TILEDB_READ, member.uri());

    throw TileDBSOMAError("Saw invalid SOMA object.");
}

bool SOMACollection::has(const std::string& key) {
    return group_->has_member(key);
}

uint64_t SOMACollection::count() const {
    return group_->get_length();
}

void SOMACollection::del(const std::string& key) {
    group_->remove_member(key);
}

std::map<std::string, std::string> SOMACollection::member_to_uri_mapping()
    const {
    return group_->member_to_uri_mapping();
};

std::unique_ptr<SOMACollection> SOMACollection::add_new_collection(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx) {
    auto member = SOMACollection::create(ctx, uri);
    group_->add_member(std::string(uri), relative, std::string(key));
    return member;
}

std::unique_ptr<SOMAExperiment> SOMACollection::add_new_experiment(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMAExperiment::create(ctx, uri, schema);
    group_->add_member(std::string(uri), relative, std::string(key));
    return member;
}

std::unique_ptr<SOMAMeasurement> SOMACollection::add_new_measurement(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMAMeasurement::create(ctx, uri, schema);
    group_->add_member(std::string(uri), relative, std::string(key));
    return member;
}

std::unique_ptr<SOMADataFrame> SOMACollection::add_new_dataframe(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMADataFrame::create(ctx, uri, schema);
    group_->add_member(std::string(uri), relative, std::string(key));
    return member;
}

std::unique_ptr<SOMADenseNDArray> SOMACollection::add_new_dense_ndarray(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMADenseNDArray::create(ctx, uri, schema);
    group_->add_member(std::string(uri), relative, std::string(key));
    return member;
}

std::unique_ptr<SOMASparseNDArray> SOMACollection::add_new_sparse_ndarray(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMASparseNDArray::create(ctx, uri, schema);
    group_->add_member(std::string(uri), relative, std::string(key));
    return member;
}

}  // namespace tiledbsoma
