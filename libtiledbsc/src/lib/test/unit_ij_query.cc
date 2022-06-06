#include <memory>
#include <stdexcept>
#include <vector>

#include <tiledbsc/ij_query.h>
#include <tiledb/tiledb>

#include <catch2/catch_test_macros.hpp>

using namespace tiledb;
using namespace tiledbsc;

// this preprocessor constant is set by CMake to *root* of source tree
const std::string src_path = TILEDBSC_SOURCE_ROOT;
const std::string data_path = src_path + "/data";

TEST_CASE("Basic IJQuery functionality") {
    Context ctx;

    std::string ref_uri = data_path + "/simple/ij_test1/ref";
    std::string tgt_uri = data_path + "/simple/ij_test1/tgt";

    auto array1 = std::make_shared<Array>(ctx, ref_uri, TILEDB_READ);
    auto array2 = std::make_shared<Array>(ctx, tgt_uri, TILEDB_READ);

    IJQuery q(array1, {"z_name", "z_idx"}, array2, "z_idx");

    // auto result = q.select_from_points<int32_t>({1,3,5,6});
};
