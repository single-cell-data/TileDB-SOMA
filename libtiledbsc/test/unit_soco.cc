#include <catch2/catch_test_macros.hpp>
#include <tiledb/tiledb>
#include <tiledbsc/tiledbsc>

#ifndef TILEDBSC_SOURCE_ROOT
#define TILEDBSC_SOURCE_ROOT "not_defined"
#endif

static const std::string root = TILEDBSC_SOURCE_ROOT;
static const std::string soco_uri = root + "/test/soco";

using namespace tiledb;
using namespace tiledbsc;

TEST_CASE("SOCO: Open arrays") {
    Config config;
    config["config.logging_level"] = "5";

    auto soco = SOMACollection::open(soco_uri, config);
    auto soma_uris = soco->list_somas();
    REQUIRE(soma_uris.size() == 2);

    for (const auto& [name, uri] : soma_uris) {
        (void)name;
        auto soma = SOMA::open(uri);
    }
}
