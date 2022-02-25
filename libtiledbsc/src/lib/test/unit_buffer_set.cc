#include <tiledb/tiledb>
#include <tiledbsc/buffer_set.h>

#include <vector>
#include <memory>

#include <catch2/catch_test_macros.hpp>

using namespace tiledb;
using namespace tiledbsc;

TEST_CASE("Basic BufferSet functionality") {
    auto b1 = tiledbsc::BufferSet(
        "foo",
        1023,
        4,
        false,
        false
    );

    REQUIRE(b1.name() == "foo");
    REQUIRE(b1.isvar() == false);
    REQUIRE(b1.isnullable() == false);
    REQUIRE(b1.data.size() == 1023 * 4);
    REQUIRE(b1.offsets == std::nullopt);
    REQUIRE(b1.validity == std::nullopt);
};

TEST_CASE("BufferSet::from_attribute") {
    auto ctx = tiledb::Context();

    auto attr = tiledb::Attribute::create<int32_t>(ctx, "a1");
    attr.set_nullable(true);
    attr.set_cell_val_num(TILEDB_VAR_NUM);

    auto buffers = BufferSet::from_attribute(attr, 21);

    REQUIRE(buffers->name() == "a1");
    REQUIRE(buffers->isvar() == true);
    REQUIRE(buffers->isnullable() == true);
    REQUIRE(buffers->data.size() == 21 * sizeof(int32_t));
    REQUIRE(buffers->offsets.value().size() == 21 * sizeof(uint64_t));
    REQUIRE(buffers->validity.value().size() == 21);
};

TEST_CASE("BufferSet::from_dimension") {
    auto ctx = tiledb::Context();
    auto dim = Dimension::create(ctx, "d1", TILEDB_STRING_ASCII, nullptr, nullptr);
    dim.set_cell_val_num(TILEDB_VAR_NUM);

    auto buffers = BufferSet::from_dimension(dim, 9);

    REQUIRE(buffers->name() == "d1");
    REQUIRE(buffers->isvar() == true);
    REQUIRE(buffers->isnullable() == false);
    REQUIRE(buffers->data.size() == 9 * sizeof(char));
    REQUIRE(buffers->offsets.value().size() == 9 * sizeof(uint64_t));
    REQUIRE(buffers->validity == std::nullopt);

    {
        size_t new_size = 100;
        buffers->resize(new_size);
        REQUIRE(buffers->data.size() == new_size);
        REQUIRE(buffers->offsets.value().size() == new_size * sizeof(uint64_t));
        REQUIRE(buffers->validity == std::nullopt);
    }
};