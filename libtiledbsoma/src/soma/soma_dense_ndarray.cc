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

#include "../cpp_api/soma_dense_ndarray.h"
#include "soma_array.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMADenseNDArray> SOMADenseNDArray::open(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::string_view name,
    std::map<std::string, std::string> platform_config,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMADenseNDArray>(
        mode,
        uri,
        name,
        std::make_shared<Context>(Config(platform_config)),
        column_names,
        batch_size,
        result_order,
        timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMADenseNDArray::SOMADenseNDArray(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::string_view name,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    array_ = std::make_unique<SOMAArray>(
        mode,
        uri,
        name,
        ctx,
        column_names,
        batch_size,
        result_order,
        timestamp);
};

void SOMADenseNDArray::close() {
    array_.get()->close();
}

std::string SOMADenseNDArray::uri() const {
    return array_.get()->uri();
}

std::shared_ptr<ArraySchema> SOMADenseNDArray::schema() const {
    return array_.get()->schema();
}

std::vector<int64_t> SOMADenseNDArray::shape() const {
    return array_.get()->shape();
}

int64_t SOMADenseNDArray::ndim() const {
    return array_.get()->ndim();
}

std::optional<std::shared_ptr<ArrayBuffers>> SOMADenseNDArray::read_next() {
    return array_.get()->read_next();
}

void SOMADenseNDArray::write(std::shared_ptr<ArrayBuffers> buffers) {
    array_.get()->write(buffers);
}

}  // namespace tiledbsoma
