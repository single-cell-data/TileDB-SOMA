/**
 * @file   unit_column_buffer.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
