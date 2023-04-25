/**
 * @file   soma_dataframe.cc
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
 *   This file defines the SOMADataFrame class.
 */

#include "../cpp_api/soma_dataframe.h"
#include <tiledb/tiledb>
#include "../utils/array_buffers.h"
#include "soma_array.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::shared_ptr<SOMADataFrame> SOMADataFrame::open(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::string_view name,
    std::map<std::string, std::string> platform_config,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_shared<SOMADataFrame>(
        mode,
        uri,
        name,
        std::make_shared<Context>(Config(platform_config)),
        column_names,
        batch_size,
        result_order,
        timestamp);
}

std::shared_ptr<SOMADataFrame> SOMADataFrame::open(
    tiledb_query_type_t mode,
    std::shared_ptr<Context> ctx,
    std::string_view uri,
    std::string_view name,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    return std::make_shared<SOMADataFrame>(
        mode,
        uri,
        name,
        ctx,
        column_names,
        batch_size,
        result_order,
        timestamp);
}
//===================================================================
//= public non-static
//===================================================================

SOMADataFrame::SOMADataFrame(
    tiledb_query_type_t mode,
    std::string_view uri,
    std::string_view name,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    array_ = std::make_shared<SOMAArray>(
        mode,
        uri,
        name,
        ctx,
        column_names,
        batch_size,
        result_order,
        timestamp);
    array_.get()->set_metadata(
        "soma_object_type", TILEDB_STRING_UTF8, 1, "SOMADataFrame");
}

void SOMADataFrame::close() {
    array_.get()->close();
}

std::string SOMADataFrame::uri() const {
    return array_.get()->uri();
}

std::shared_ptr<Context> SOMADataFrame::ctx() {
    return array_.get()->ctx();
}

std::shared_ptr<ArraySchema> SOMADataFrame::schema() const {
    return array_.get()->schema();
}

std::vector<std::string> SOMADataFrame::index_column_names() const {
    return array_.get()->dimension_names();
}

int64_t SOMADataFrame::count() const {
    return array_.get()->ndim();
}

std::optional<std::shared_ptr<ArrayBuffers>> SOMADataFrame::read_next() {
    return array_.get()->read_next();
}

void SOMADataFrame::write(std::shared_ptr<ArrayBuffers> buffers) {
    array_.get()->write(buffers);
}

}  // namespace tiledbsoma
