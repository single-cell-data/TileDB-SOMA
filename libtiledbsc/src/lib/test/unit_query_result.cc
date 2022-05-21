#include <filesystem>
#include <memory>
#include <stdexcept>
#include <vector>

#include <tiledbsc/query_result.h>
#include <tiledbsc/sc_arrowio.h>
#include <tiledbsc/util.h>
#include <tiledb/tiledb>

using namespace std;
using namespace tiledbsc;

#include <catch2/catch_test_macros.hpp>

using namespace tiledbsc;

namespace {
using namespace tiledbsc;

ResultBuffers create_data1() {
    std::vector<string> data{
        "", "abcd", "ef", "ghijk", "lmno", "p", "", "q", "rstu", "vwxyz", ""};

    std::vector<uint8_t> validity_vec_tmp{0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1};
    std::vector<byte> validity_vec;
    validity_vec.resize(11);
    std::transform(
        validity_vec_tmp.begin(),
        validity_vec_tmp.end(),
        validity_vec.begin(),
        [](uint8_t v) -> byte { return (byte)v; });

    auto [data_buf, offsets_buf] = util::to_varlen_buffers(data);

    auto bufferset_b = BufferSet::from_data(
        "b", TILEDB_STRING_ASCII, 1, data_buf, offsets_buf, validity_vec);

    return ResultBuffers{
        std::make_pair(std::string("b"), std::move(bufferset_b))};
}

};  // anonymous namespace

TEST_CASE("QueryResult constructor") {
    auto result_buffers = create_data1();
    auto query_result = QueryResult{std::move(result_buffers)};

    auto arrow_array = query_result.to_arrow("b");
};
