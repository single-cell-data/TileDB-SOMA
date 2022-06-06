#include <bitset>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <tiledbsc/util.h>

#include <catch2/catch_test_macros.hpp>

using namespace std;
using namespace tiledbsc;

// work-around for difference between c++2a and c++20
//   span has index_type in 2a which doesn't match size_t in 20
namespace {
using SPAN_T = std::span<const byte>;
template <typename T, typename = void>
struct check_index_type {
    using type = size_t;
};

template <typename T>
struct check_index_type<T, void_t<typename T::index_type>> {
    using type = typename T::index_type;
};

// ditto for bar, baz, etc.

template <typename T>
struct span_index_type : check_index_type<T> {};

};  // namespace

TEST_CASE("Test to_varlen_buffers") {
    vector<string> data{
        "", "abcd", "ef", "ghijk", "lmno", "p", "", "q", "rstu", "vwxyz", ""};

    size_t data_size = 0;
    std::vector<byte> data_raw;

    std::for_each(
        data.begin(), data.end(), [&data_size, &data_raw](std::string& s) {
            data_size += s.size();
            const span<const byte> s_span{
                reinterpret_cast<const byte*>(s.data()),
                static_cast<span_index_type<SPAN_T>::type>(s.size())};
            data_raw.insert(data_raw.end(), s_span.begin(), s_span.end());
        });

    auto [data_buf, offsets_buf] = util::to_varlen_buffers(data);

    REQUIRE(offsets_buf.size() == data.size() + 1);
    REQUIRE(!memcmp(data_buf.data(), data_raw.data(), data_buf.size()));
    REQUIRE(
        offsets_buf ==
        vector<uint64_t>{0, 0, 4, 6, 11, 15, 16, 16, 17, 21, 26, 26});
}

TEST_CASE("Util: Arrow bytemap to bitmap conversion", "[util][arrow][bitmap]") {
    {
        std::vector<uint8_t> validity{1};

        // REQUIRE(util::bytemap_to_bitmap_inplace(validity) == 0);
    }

    {
        std::vector<uint8_t> validity{1, 0, 1, 0, 1, 0, 0, 0};

        REQUIRE(util::bytemap_to_bitmap_inplace(validity) == 5);

        REQUIRE(std::bitset<8>(validity[0]) == 0b00010101);
    }

    {
        std::vector<uint8_t> validity{1,
                                      0,
                                      1,
                                      0,
                                      1,
                                      0,
                                      0,
                                      1,  // break
                                      0,
                                      1,
                                      1};

        REQUIRE(util::bytemap_to_bitmap_inplace(validity) == 5);

        REQUIRE(std::bitset<8>(validity[0]) == 0b10010101);
        REQUIRE(std::bitset<8>(validity[1]) == 0b00000110);
    }

    {
        std::vector<uint8_t> validity{0,
                                      0,
                                      0,
                                      1,
                                      0,
                                      1,
                                      1,
                                      1,  // break
                                      0,
                                      1,
                                      1,
                                      0,
                                      0,
                                      0,
                                      1,
                                      0,  // break
                                      1};

        REQUIRE(util::bytemap_to_bitmap_inplace(validity) == 9);

        REQUIRE(std::bitset<8>(validity[0]) == 0b11101000);
        REQUIRE(std::bitset<8>(validity[1]) == 0b01000110);
        REQUIRE(std::bitset<8>(validity[2]) == 0b00000001);
    }
}