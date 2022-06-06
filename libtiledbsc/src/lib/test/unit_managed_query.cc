#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <tiledbsc/managed_query.h>
#include <tiledb/tiledb>

#include <catch2/catch_test_macros.hpp>

// this preprocessor constant is set by CMake to *root* of source tree
const std::string src_path = TILEDBSC_SOURCE_ROOT;

namespace {

const std::unordered_map<std::string, std::string> global_conf{
    {"sm.query.dense.reader", "legacy"}  // SC-17113
};

void init_conf(tiledb::Config& conf) {
    for (auto& c : global_conf) {
        conf[c.first] = c.second;
    }
}

tiledb::Context make_ctx() {
    tiledb::Config config;
    init_conf(config);
    return tiledb::Context(config);
}

};  // end anonymous namespace

TEST_CASE("Basic test of ManagedQuery execution and results") {
    // path to data within the *source* tree
    auto data_path = src_path + "/data/";
    auto array_path = data_path + "/simple/dim-uint64_attr-int64_10cells/";

    auto ctx = make_ctx();

    auto array = std::make_shared<tiledb::Array>(ctx, array_path, TILEDB_READ);

    auto mq = tiledbsc::ManagedQuery(array);

    mq.select_points(0, std::vector<uint64_t>{1});

    auto result = mq.execute();

    auto data = result->get("");
    REQUIRE(data.data_.size() == 8);  // bytes
    auto dim = result->get("__dim_0");
    REQUIRE(dim.data_.size() == 8);
};

TEST_CASE("ManagedQuery string attribute test") {
    // path to data within the *source* tree
    auto data_path = src_path + "/data/";
    auto array_path = data_path + "/simple/dim-uint64_attr-str_26cells/";

    auto ctx = make_ctx();

    auto array = std::make_shared<tiledb::Array>(ctx, array_path, TILEDB_READ);

    auto mq = tiledbsc::ManagedQuery(array);
    mq.select_points(0, std::vector<uint64_t>{1, 3, 6});
    auto result = mq.execute();

    REQUIRE_THROWS(result->get("foobar"));

    auto& buf = result->get("");
    REQUIRE(buf.data_.size() == 13);  // char
    REQUIRE(std::string_view((char*)buf.data_.data(), 13) == "bbddddggggggg");

    auto dim = result->get("__dim_0");
    REQUIRE(dim.data_.size() == 24);  // uint64

    REQUIRE(result->names() == std::vector<std::string>({"", "__dim_0"}));
    REQUIRE(result->nbuffers() == 2);
};