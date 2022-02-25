#include <tiledb/tiledb>
#include <tiledbsc/managed_query.h>

#include <memory>
#include <filesystem>
#include <vector>

#include <catch2/catch_test_macros.hpp>

// this preprocessor constant is set by CMake to *root* of source tree
const std::string src_path = TILEDBSC_SOURCE_ROOT;

TEST_CASE("Basic test of ManagedQuery execution and results") {
    // path to data within the *source* tree
    auto data_path = src_path + "/data/";
    auto array_path = data_path + "/simple/dim-uint64_attr-int64_10cells/";

    auto ctx = tiledb::Context();
    auto array = std::make_shared<tiledb::Array>(ctx, array_path, TILEDB_READ);

    auto mq = tiledbsc::ManagedQuery(array);

    mq.select_points(0, std::vector<uint64_t>{1});

    auto result = mq.execute();

    auto data = result->get("");
    REQUIRE(data != std::nullopt);
    REQUIRE(data.value()->data.size() == 8); // bytes
    auto dim = result->get("__dim_0");
    REQUIRE(dim != std::nullopt);
    REQUIRE(dim.value()->data.size() == 8);
};

TEST_CASE("ManagedQuery string attribute test") {
    // path to data within the *source* tree
    auto data_path = src_path + "/data/";
    auto array_path = data_path + "/simple/dim-uint64_attr-str_26cells/";

    auto ctx = tiledb::Context();
    auto array = std::make_shared<tiledb::Array>(ctx, array_path, TILEDB_READ);

    auto mq = tiledbsc::ManagedQuery(array);
    mq.select_points(0, std::vector<uint64_t>{1,3,6});
    auto result = mq.execute();

    auto data = result->get("");
    REQUIRE(data != std::nullopt);
    REQUIRE(data.value()->data.size() == 13); // char
    auto dim = result->get("__dim_0");
    REQUIRE(dim != std::nullopt);
    REQUIRE(dim.value()->data.size() == 24); // uint64

    REQUIRE(std::string((char*)data.value()->data.data(), 13) == "bbddddggggggg");

    REQUIRE(result->names() == std::vector<std::string>({"", "__dim_0"}));
    REQUIRE(result->nbuffers() == 2);
};