#include <tiledbsc/buffer_set.h>
#include <tiledb/tiledb>

#include <memory>
#include <stdexcept>  // for windows error C2039
#include <vector>

#include <catch2/catch_test_macros.hpp>

using namespace tiledb;

const std::string src_path = TILEDBSC_SOURCE_ROOT;

TEST_CASE("placeholder") {}