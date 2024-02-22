/**
 * @file   soma_sparse_ndarray.cc
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
 *   This file defines the SOMASparseNDArray class.
 */

#include "soma_sparse_ndarray.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMASparseNDArray> SOMASparseNDArray::create(
    std::string_view uri,
    ArraySchema schema,
    std::map<std::string, std::string> platform_config) {
    return SOMASparseNDArray::create(
        uri, schema, std::make_shared<Context>(Config(platform_config)));
}

std::unique_ptr<SOMASparseNDArray> SOMASparseNDArray::create(
    std::string_view uri, ArraySchema schema, std::shared_ptr<Context> ctx) {
    if (schema.array_type() != TILEDB_SPARSE)
        throw TileDBSOMAError("ArraySchema must be set to sparse.");

    SOMAArray::create(ctx, uri, schema, "SOMASparseNDArray");
    return SOMASparseNDArray::open(uri, OpenMode::read, ctx);
}

std::unique_ptr<SOMASparseNDArray> SOMASparseNDArray::open(
    std::string_view uri,
    OpenMode mode,
    std::map<std::string, std::string> platform_config,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return SOMASparseNDArray::open(
        uri,
        mode,
        std::make_shared<Context>(Config(platform_config)),
        column_names,
        result_order,
        timestamp);
}

std::unique_ptr<SOMASparseNDArray> SOMASparseNDArray::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_unique<SOMASparseNDArray>(
        mode, uri, ctx, column_names, result_order, timestamp);
}

bool SOMASparseNDArray::exists(std::string_view uri) {
    try {
        return "SOMASparseNDArray" ==
               SOMAObject::open(uri, OpenMode::read)->type();
    } catch (TileDBSOMAError& e) {
        return false;
    }
}

//===================================================================
//= public non-static
//===================================================================

std::unique_ptr<ArrowSchema> SOMASparseNDArray::schema() const {
    return this->arrow_schema();
}
}  // namespace tiledbsoma
