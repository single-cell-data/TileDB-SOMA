/**
 * @file   unit_column_buffer.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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
 * This file manages unit tests for column buffers
 */

#include <catch2/catch_test_macros.hpp>
#include <tiledb/tiledb>
#include <tiledbsoma/tiledbsoma>

using namespace tiledb;
using namespace tiledbsoma;

#ifndef TILEDBSOMA_SOURCE_ROOT
#define TILEDBSOMA_SOURCE_ROOT "not_defined"
#endif

const std::string src_path = TILEDBSOMA_SOURCE_ROOT;

namespace {

/**
 * @brief Create an array and return array opened in read mode.
 *
 * @param uri Array uri
 * @param ctx TileDB context
 * @return std::shared_ptr<Array>
 */
static std::shared_ptr<Array> create_array(
    const std::string& uri, Context& ctx) {
    // delete array if it exists
    auto vfs = VFS(ctx);
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    ArraySchema schema(ctx, TILEDB_SPARSE);
    auto dim = Dimension::create(
        ctx, "d1", TILEDB_STRING_ASCII, nullptr, nullptr);
    dim.set_cell_val_num(TILEDB_VAR_NUM);

    Domain domain(ctx);
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int32_t>(ctx, "a1");
    attr.set_nullable(true);
    attr.set_cell_val_num(TILEDB_VAR_NUM);
    schema.add_attribute(attr);

    Array::create(uri, std::move(schema));
    return std::make_shared<Array>(ctx, uri, TILEDB_READ);
}

};  // namespace

TEST_CASE("ColumnBuffer: Create from array") {
    std::string uri = "mem://unit-test-array";
    auto ctx = Context();
    auto array = create_array(uri, ctx);

    {
        auto buffers = ColumnBuffer::create(array, "d1");
        REQUIRE(buffers->name() == "d1");
        REQUIRE(buffers->is_var() == true);
        REQUIRE(buffers->is_nullable() == false);
    }

    {
        auto buffers = ColumnBuffer::create(array, "a1");
        REQUIRE(buffers->name() == "a1");
        REQUIRE(buffers->is_var() == true);
        REQUIRE(buffers->is_nullable() == true);
    }
}
