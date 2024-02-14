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
    std::string_view uri, std::map<std::string, std::string> platform_config) {
    return SOMACollection::create(
        uri, std::make_shared<Context>(Config(platform_config)));
}

std::unique_ptr<SOMACollection> SOMACollection::create(
    std::string_view uri, std::shared_ptr<Context> ctx) {
    SOMAGroup::create(ctx, uri, "SOMACollection");
    return SOMACollection::open(uri, OpenMode::read, ctx);
}

std::unique_ptr<SOMACollection> SOMACollection::open(
    std::string_view uri,
    OpenMode mode,
    std::map<std::string, std::string> platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return SOMACollection::open(
        uri,
        mode,
        std::make_shared<Context>(Config(platform_config)),
        timestamp);
}

std::unique_ptr<SOMACollection> SOMACollection::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<Context> ctx,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMACollection>(mode, uri, ctx, timestamp);
}

//===================================================================
//= public non-static
//===================================================================

void SOMACollection::set(
    std::string_view uri, URIType uri_type, const std::string& key) {
    this->add_member(std::string(uri), uri_type, key);
}

std::shared_ptr<SOMAObject> SOMACollection::get(const std::string& key) {
    auto member = this->get_member(key);
    std::optional<std::string> soma_object_type = this->type();

    if (!soma_object_type)
        throw TileDBSOMAError("Saw invalid SOMA object.");

    if (soma_object_type->compare("SOMACollection") == 0)
        return SOMACollection::open(member.uri(), OpenMode::read);
    else if (soma_object_type->compare("SOMAExperiment") == 0)
        return SOMAExperiment::open(member.uri(), OpenMode::read);
    else if (soma_object_type->compare("SOMAMeasurement") == 0)
        return SOMAMeasurement::open(member.uri(), OpenMode::read);
    else if (soma_object_type->compare("SOMADataFrame") == 0)
        return SOMADataFrame::open(member.uri(), OpenMode::read);
    else if (soma_object_type->compare("SOMASparseNDArray") == 0)
        return SOMASparseNDArray::open(member.uri(), OpenMode::read);
    else if (soma_object_type->compare("SOMADenseNDArray") == 0)
        return SOMADenseNDArray::open(member.uri(), OpenMode::read);

    throw TileDBSOMAError("Saw invalid SOMA object.");
}

bool SOMACollection::has(const std::string& key) {
    return this->has_member(key);
}

uint64_t SOMACollection::count() const {
    return this->get_length();
}

void SOMACollection::del(const std::string& key) {
    this->remove_member(key);
}

std::unique_ptr<SOMACollection> SOMACollection::add_new_collection(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<Context> ctx) {
    auto member = SOMACollection::create(uri, ctx);
    this->add_member(std::string(uri), uri_type, std::string(key));
    children_[std::string(key)] = std::make_shared<SOMACollection>(*member);
    return member;
}

std::unique_ptr<SOMAExperiment> SOMACollection::add_new_experiment(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMAExperiment::create(uri, schema, ctx);
    this->add_member(std::string(uri), uri_type, std::string(key));
    children_[std::string(key)] = std::make_shared<SOMAExperiment>(*member);
    return member;
}

std::unique_ptr<SOMAMeasurement> SOMACollection::add_new_measurement(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMAMeasurement::create(uri, schema, ctx);
    this->add_member(std::string(uri), uri_type, std::string(key));
    children_[std::string(key)] = std::make_shared<SOMAMeasurement>(*member);
    return member;
}

std::unique_ptr<SOMADataFrame> SOMACollection::add_new_dataframe(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMADataFrame::create(uri, schema, ctx);
    this->add_member(std::string(uri), uri_type, std::string(key));
    children_[std::string(key)] = std::make_shared<SOMADataFrame>(*member);
    return member;
}

std::unique_ptr<SOMADenseNDArray> SOMACollection::add_new_dense_ndarray(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMADenseNDArray::create(uri, schema, ctx);
    this->add_member(std::string(uri), uri_type, std::string(key));
    children_[std::string(key)] = std::make_shared<SOMADenseNDArray>(*member);
    return member;
}

std::unique_ptr<SOMASparseNDArray> SOMACollection::add_new_sparse_ndarray(
    std::string_view key,
    std::string_view uri,
    URIType uri_type,
    std::shared_ptr<Context> ctx,
    ArraySchema schema) {
    auto member = SOMASparseNDArray::create(uri, schema, ctx);
    this->add_member(std::string(uri), uri_type, std::string(key));
    children_[std::string(key)] = std::make_shared<SOMASparseNDArray>(*member);
    return member;
}

}  // namespace tiledbsoma
