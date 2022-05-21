#include <memory>
#include <stdexcept>  // for windows
#include <vector>

#include <tiledbsc/buffer_set.h>
#include <tiledbsc/util.h>
#include <tiledb/tiledb>

#include <catch2/catch_test_macros.hpp>

using namespace tiledb;
using namespace tiledbsc;

TEST_CASE("BufferSet::from_nelem") {
    auto b1 = tiledbsc::BufferSet::alloc(
        "foo", TILEDB_INT32, 1023, 4, false, false);

    REQUIRE(b1.name() == "foo");
    REQUIRE(b1.isvar() == false);
    REQUIRE(b1.isnullable() == false);
    REQUIRE(b1.data<byte>().size() == 1023 * 4);
    REQUIRE(b1.data_.size() == 1023 * 4);
    REQUIRE(b1.offsets_ == std::nullopt);
    REQUIRE(b1.validity_ == std::nullopt);
};

TEST_CASE("BufferSet::from_attribute") {
    auto ctx = tiledb::Context();

    auto attr = tiledb::Attribute::create<int32_t>(ctx, "a1");
    attr.set_nullable(true);
    attr.set_cell_val_num(TILEDB_VAR_NUM);

    auto buffers = BufferSet::from_attribute(attr, 21);

    REQUIRE(buffers.name() == "a1");
    REQUIRE(buffers.isvar() == true);
    REQUIRE(buffers.isnullable() == true);
    REQUIRE(buffers.data_.size() == 21 * sizeof(int32_t));
    REQUIRE(buffers.offsets_.value().size() == 22);
    REQUIRE(buffers.validity_.value().size() == 21);
};

TEST_CASE("BufferSet::from_dimension") {
    auto ctx = tiledb::Context();
    auto dim = Dimension::create(
        ctx, "d1", TILEDB_STRING_ASCII, nullptr, nullptr);
    dim.set_cell_val_num(TILEDB_VAR_NUM);

    auto buffers = BufferSet::from_dimension(dim, 9);

    REQUIRE(buffers.name() == "d1");
    REQUIRE(buffers.isvar() == true);
    REQUIRE(buffers.isnullable() == false);
    REQUIRE(buffers.data_.size() == 9 * sizeof(char));
    REQUIRE(buffers.offsets_.value().size() == 10);
    REQUIRE(buffers.validity_ == std::nullopt);

    {
        size_t new_size = 100;
        buffers.resize(new_size);
        REQUIRE(buffers.data_.size() == new_size);
        REQUIRE(buffers.offsets_.value().size() == 10);
        REQUIRE(buffers.validity_ == std::nullopt);

        buffers.resize(new_size, new_size);
        REQUIRE(buffers.data_.size() == new_size);
        REQUIRE(buffers.offsets_.value().size() == 101);
        REQUIRE(buffers.validity_ == std::nullopt);
    }
};

TEST_CASE("BufferSet::from_data") {
    std::vector<string> data_orig{
        "", "abcd", "ef", "ghijk", "lmno", "p", "", "q", "rstu", "vwxyz", ""};

    auto&& [data_buf, offsets_buf] = util::to_varlen_buffers(data_orig);

    auto b = BufferSet::from_data(
        "b", TILEDB_STRING_ASCII, 1, data_buf, offsets_buf);

    CHECK_THROWS(b.validity());
    REQUIRE(b.data<byte>().size() == 26);
    REQUIRE(b.offsets().size() == 12);
    REQUIRE(b.offsets()[0] == 0);
    REQUIRE(b.offsets()[1] == 0);
    REQUIRE(b.offsets()[2] == 4);
    REQUIRE(b.offsets()[10] == 26);

    std::vector<uint8_t> validity_vec_tmp{0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1};
    std::vector<byte> validity_vec;
    validity_vec.resize(11);
    std::transform(
        validity_vec_tmp.begin(),
        validity_vec_tmp.end(),
        validity_vec.begin(),
        [](uint8_t v) -> byte { return (byte)v; });
    {
        auto b2 = BufferSet::from_data(
            "b", TILEDB_STRING_ASCII, 1, data_buf, offsets_buf, validity_vec);

        REQUIRE(b2.validity().size() == 11);
    }
}
