#include <catch2/catch_test_macros.hpp>
#include <tiledb/tiledb>
#include <tiledbsc/tiledbsc>

#define VERBOSE 0

#ifndef TILEDBSC_SOURCE_ROOT
#define TILEDBSC_SOURCE_ROOT "not_defined"
#endif

static const std::string root = TILEDBSC_SOURCE_ROOT;
static const std::string soma_uri = root + "/test/soco/pbmc3k_processed";

using namespace tiledb;
using namespace tiledbsc;

int soma_num_cells(MultiArrayBuffers& soma) {
    return soma.begin()->second.begin()->second->size();
}

TEST_CASE("SOMA: Open arrays") {
    if (VERBOSE) {
        LOG_CONFIG("debug");
    }

    Config config;
    // config.logging_level"] = "5";

    auto soma = SOMA::open(soma_uri, config);
    auto array_uris = soma->list_arrays();
    REQUIRE(array_uris.size() == 19);

    for (const auto& [name, uri] : array_uris) {
        (void)uri;
        auto array = soma->open_array(name);
    }
}

TEST_CASE("SOMA: Full query") {
    auto soma = SOMA::open(soma_uri);
    auto sq = soma->query();

    size_t total_cells = 0;
    while (auto results = sq->next_results()) {
        auto num_cells = soma_num_cells(*results);
        total_cells += num_cells;
    }
    REQUIRE(total_cells == 4848644);
}

TEST_CASE("SOMA: Sliced query (obs)") {
    auto soma = SOMA::open(soma_uri);
    auto sq = soma->query();
    auto ctx = soma->context();

    // Set obs query condition
    std::string obs_attr = "louvain";
    std::string obs_val = "B cells";
    auto obs_qc = QueryCondition::create(*ctx, obs_attr, obs_val, TILEDB_EQ);
    std::vector<std::string> obs_cols = {obs_attr};
    sq->set_obs_condition(obs_qc);
    sq->select_obs_attrs(obs_cols);

    size_t total_cells = 0;
    while (auto results = sq->next_results()) {
        auto num_cells = soma_num_cells(*results);
        total_cells += num_cells;
    }
    REQUIRE(total_cells == 628596);
}

TEST_CASE("SOMA: Sliced query (var)") {
    auto soma = SOMA::open(soma_uri);
    auto sq = soma->query();
    auto ctx = soma->context();

    // Set var query condition
    std::string var_attr = "n_cells";
    uint64_t var_val = 50;
    auto var_qc = QueryCondition::create<uint64_t>(
        *ctx, var_attr, var_val, TILEDB_LT);
    std::vector<std::string> var_cols = {var_attr};
    sq->set_var_condition(var_qc);
    sq->select_var_attrs(var_cols);

    size_t total_cells = 0;
    while (auto results = sq->next_results()) {
        auto num_cells = soma_num_cells(*results);
        total_cells += num_cells;
    }
    REQUIRE(total_cells == 1308448);
}

TEST_CASE("SOMA: Sliced query (select ids)") {
    auto soma = SOMA::open(soma_uri);
    auto sq = soma->query();

    std::vector<std::string> obs_ids = {
        "AAACATACAACCAC-1", "AAACATTGATCAGC-1", "TTTGCATGCCTCAC-1"};
    std::vector<std::string> var_ids = {"AAGAB", "AAR2", "ZRANB3"};
    sq->select_obs_ids(obs_ids);
    sq->select_var_ids(var_ids);

    size_t total_cells = 0;
    while (auto results = sq->next_results()) {
        auto num_cells = soma_num_cells(*results);
        total_cells += num_cells;
    }
    REQUIRE(total_cells == 9);
}
