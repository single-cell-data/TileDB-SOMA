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
#include "../cpp_api/soma_dataframe.h"
#include "../cpp_api/soma_dense_ndarray.h"
#include "../cpp_api/soma_experiment.h"
#include "../cpp_api/soma_measurement.h"
#include "../cpp_api/soma_sparse_ndarray.h"
#include "soma_array.h"
#include "soma_group.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::shared_ptr<SOMACollection> SOMACollection::open(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::map<std::string, std::string> platform_config) {
    return std::make_shared<SOMACollection>(
        mode, uri, std::make_shared<Context>(Config(platform_config)));
}

std::shared_ptr<SOMACollection> SOMACollection::open(
    tiledb_query_type_t mode,
    std::shared_ptr<Context> ctx,
    std::string_view uri) {
    return std::make_shared<SOMACollection>(mode, uri, ctx);
}

//===================================================================
//= public non-static
//===================================================================

SOMACollection::SOMACollection(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::shared_ptr<Context> ctx)
    : ctx_(ctx) {
    group_ = std::make_shared<SOMAGroup>(mode, uri, ctx_);
    group_.get()->set_metadata(
        "soma_object_type", TILEDB_STRING_UTF8, 1, "SOMACollection");
}

void SOMACollection::close() {
    group_.get()->close();
}

std::string SOMACollection::type() const {
    tiledb_datatype_t value_type;
    uint32_t value_num;
    const void* value;

    group_.get()->get_metadata(
        "soma_object_type", &value_type, &value_num, &value);

    return std::string(static_cast<const char*>(value), value_num);
}

std::string SOMACollection::uri() const {
    return group_.get()->uri();
}

std::shared_ptr<Context> SOMACollection::ctx() {
    return group_.get()->ctx();
}

void SOMACollection::set(const std::string& key, SOMAObject& object) {
    // TODO need to check if URI of object is relative
    group_.get()->add_member(object.uri(), false, key);
}

std::shared_ptr<SOMAObject> SOMACollection::get(const std::string& key) {
    auto member = group_.get()->get_member(key);
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
    return group_.get()->has_member(key);
}

uint64_t SOMACollection::count() const {
    return group_.get()->get_length();
}

void SOMACollection::del(const std::string& key) {
    group_.get()->remove_member(key);
}

std::map<std::string, std::string> SOMACollection::member_to_uri_mapping()
    const {
    return group_.get()->member_to_uri_mapping();
};

SOMACollection SOMACollection::add_new_collection(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx) {
    auto member = SOMACollection(TILEDB_READ, uri, ctx);
    group_.get()->add_member(std::string(uri), relative, std::string(key));
    return member;
}

// SOMAExperiment SOMACollection::add_new_experiment(
//     std::string_view key,
//     std::string_view uri,
//     bool relative,
//     std::shared_ptr<Context> ctx,
//     SOMADataFrame& obs,
//     SOMACollection& ms) {
//     auto member = SOMAExperiment(TILEDB_READ, uri, ctx, obs, ms);
//     group_.get()->add_member(std::string(uri), relative, std::string(key));
//     return member;
// }

// SOMAMeasurement SOMACollection::add_new_measurement(
//     std::string_view key,
//     std::string_view uri,
//     bool relative,
//     std::shared_ptr<Context> ctx,
//     SOMADataFrame& var,
//     SOMACollection& X,
//     SOMACollection& obsm,
//     SOMACollection& obsp,
//     SOMACollection& varm,
//     SOMACollection& varp) {
//     auto member = SOMAMeasurement(
//         TILEDB_READ, uri, ctx, var, X, obsm, obsp, varm, varp);
//     group_.get()->add_member(std::string(uri), relative, std::string(key));
//     return member;
// }

SOMADataFrame SOMACollection::add_new_dataframe(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto member = SOMADataFrame(
        TILEDB_READ,
        uri,
        key,
        ctx,
        column_names,
        batch_size,
        result_order,
        timestamp);
    group_.get()->add_member(std::string(uri), relative, std::string(key));
    return member;
}

SOMADenseNDArray SOMACollection::add_new_dense_ndarray(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto member = SOMADenseNDArray(
        TILEDB_READ,
        uri,
        key,
        ctx,
        column_names,
        batch_size,
        result_order,
        timestamp);
    group_.get()->add_member(std::string(uri), relative, std::string(key));
    return member;
}

SOMASparseNDArray SOMACollection::add_new_sparse_ndarray(
    std::string_view key,
    std::string_view uri,
    bool relative,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto member = SOMASparseNDArray(
        TILEDB_READ,
        uri,
        key,
        ctx,
        column_names,
        batch_size,
        result_order,
        timestamp);
    group_.get()->add_member(std::string(uri), relative, std::string(key));
    return member;
}

}  // namespace tiledbsoma
