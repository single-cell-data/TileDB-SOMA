#include <catch2/catch_test_macros.hpp>
#include <tiledb/tiledb>
#include <tiledbsc/tiledbsc>

using namespace tiledb;
using namespace tiledbsc;

#ifndef TILEDBSC_SOURCE_ROOT
#define TILEDBSC_SOURCE_ROOT "not_defined"
#endif

const std::string src_path = TILEDBSC_SOURCE_ROOT;

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

    Array::create(uri, schema);
    return std::make_shared<Array>(ctx, uri, TILEDB_READ);
}

TEST_CASE("ColumnBuffer: Create from attribute") {
    auto ctx = Context();

    auto attr = Attribute::create<int32_t>(ctx, "a1");
    attr.set_nullable(true);
    attr.set_cell_val_num(TILEDB_VAR_NUM);

    auto buffers = ColumnBuffer::create(attr, 21);

    REQUIRE(buffers->name() == "a1");
    REQUIRE(buffers->is_var() == true);
    REQUIRE(buffers->is_nullable() == true);
    REQUIRE(buffers->data<int32_t>().size() == 21);
    REQUIRE(buffers->offsets().size() == 22);
    REQUIRE(buffers->validity().size() == 21);
}

TEST_CASE("ColumnBuffer: Create from dimension") {
    auto ctx = Context();
    auto dim = Dimension::create(
        ctx, "d1", TILEDB_STRING_ASCII, nullptr, nullptr);
    dim.set_cell_val_num(TILEDB_VAR_NUM);

    auto buffers = ColumnBuffer::create(dim, 9);

    REQUIRE(buffers->name() == "d1");
    REQUIRE(buffers->is_var() == true);
    REQUIRE(buffers->is_nullable() == false);
    REQUIRE(buffers->data<char>().size() == 9);
    REQUIRE(buffers->offsets().size() == 10);
    REQUIRE(buffers->validity().size() == 0);
}

TEST_CASE("ColumnBuffer: Create from array") {
    std::string uri = "mem://unit-test-array";
    auto ctx = Context();
    auto array = create_array(uri, ctx);

    {
        auto buffers = ColumnBuffer::create(array, "d1", 9);
        REQUIRE(buffers->name() == "d1");
        REQUIRE(buffers->is_var() == true);
        REQUIRE(buffers->is_nullable() == false);
        REQUIRE(buffers->data<char>().size() == 9);
        REQUIRE(buffers->offsets().size() == 10);
        REQUIRE(buffers->validity().size() == 0);
    }

    {
        auto buffers = ColumnBuffer::create(array, "a1", 21);
        REQUIRE(buffers->name() == "a1");
        REQUIRE(buffers->is_var() == true);
        REQUIRE(buffers->is_nullable() == true);
        REQUIRE(buffers->data<int32_t>().size() == 21);
        REQUIRE(buffers->offsets().size() == 22);
        REQUIRE(buffers->validity().size() == 21);
    }
}
