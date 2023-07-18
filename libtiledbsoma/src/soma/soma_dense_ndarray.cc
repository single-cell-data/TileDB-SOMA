/**
 * @file   soma_dense_ndarray.cc
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
 *   This file defines the SOMADenseNDArray class.
 */

#include "soma_dense_ndarray.h"
#include "soma_array.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMADenseNDArray> SOMADenseNDArray::create(
    std::shared_ptr<Context> ctx, std::string_view uri, ArraySchema schema) {
    if (schema.array_type() != TILEDB_DENSE)
        throw TileDBSOMAError("ArraySchema must be set to dense.");

    SOMAArray::create(ctx, uri, schema, "SOMADenseNDArray");
    return std::make_unique<SOMADenseNDArray>(
        TILEDB_READ, uri, ctx, std::vector<std::string>(), std::nullopt);
}

std::unique_ptr<SOMADenseNDArray> SOMADenseNDArray::open(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::vector<std::string> column_names,
    std::map<std::string, std::string> platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto ctx = std::make_shared<Context>(Config(platform_config));
    return std::make_unique<SOMADenseNDArray>(
        mode, uri, ctx, column_names, timestamp);
}

std::unique_ptr<SOMADenseNDArray> SOMADenseNDArray::open(
    tiledb_query_type_t mode,
    std::shared_ptr<Context> ctx,
    std::string_view uri,
    std::vector<std::string> column_names,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMADenseNDArray>(
        mode, uri, ctx, column_names, timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMADenseNDArray::SOMADenseNDArray(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    array_ = std::make_shared<SOMAArray>(
        mode,
        uri,
        "unnamed",  // name
        ctx,
        column_names,
        "auto",  // batch_size,
        "auto",  // result_order,
        timestamp);
    array_->reset();
    array_->submit();
}

void SOMADenseNDArray::open(
    tiledb_query_type_t mode,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    array_->open(mode, timestamp);
    array_->reset();
    array_->submit();
}

void SOMADenseNDArray::close() {
    array_->close();
}

const std::string SOMADenseNDArray::uri() const {
    return array_->uri();
}

std::shared_ptr<Context> SOMADenseNDArray::ctx() {
    return array_->ctx();
}

std::shared_ptr<ArraySchema> SOMADenseNDArray::schema() const {
    return array_->schema();
}

std::vector<int64_t> SOMADenseNDArray::shape() const {
    return array_->shape();
}

int64_t SOMADenseNDArray::ndim() const {
    return array_->ndim();
}

std::optional<std::shared_ptr<ArrayBuffers>> SOMADenseNDArray::read_next() {
    return array_->read_next();
}

void SOMADenseNDArray::write(std::shared_ptr<ArrayBuffers> buffers) {
    array_->reset();
    array_->submit();
    array_->write(buffers);
}

}  // namespace tiledbsoma
