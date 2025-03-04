/**
 * @file   cli.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file is currently a sandbox for C++ API experiments
 */

#include <format>
#include "soma/enums.h"
#include "soma/soma_array.h"
#include "utils/arrow_adapter.h"
#include "utils/logger.h"

using namespace tiledbsoma;

// [[Rcpp::export]]
void test_sdf(const std::string& uri) {
    auto ctx = std::make_shared<SOMAContext>();

    std::map<std::string, std::string> config;
    // Control buffer sizes, similar to tiledb-py
    // config["soma.init_buffer_bytes"] = "4294967296";

    // Control core memory usage
    // config["sm.mem.total_budget"] = "1118388608";

    // Read all values from the obs array
    auto obs = SOMAArray::open(OpenMode::read, uri + "/obs", ctx);
    auto obs_mq = ManagedQuery(*obs, ctx->tiledb_ctx(), "obs");
    auto obs_data = obs_mq.read_next();

    // Read all values from the var array
    auto var = SOMAArray::open(OpenMode::read, uri + "/ms/RNA/var", ctx);
    auto var_mq = ManagedQuery(*var, ctx->tiledb_ctx(), "var");
    auto var_data = var_mq.read_next();

    // Check if obs and var reads are complete
    if (obs_mq.results_complete() && var_mq.results_complete()) {
        LOG_INFO("var and obs queries are complete");
    }

    // Read all values from the X/data array
    auto x_data = SOMAArray::open(
        OpenMode::read,
        uri + "/ms/RNA/X/data",
        std::make_shared<SOMAContext>(config));
    auto x_mq = ManagedQuery(*x_data, ctx->tiledb_ctx(), "X/data");

    int batches = 0;
    int total_num_rows = 0;

    // Handle incomplete queries
    while (auto batch = x_mq.read_next()) {
        batches++;
        total_num_rows += (*batch)->num_rows();
    }

    LOG_INFO(std::format("X/data rows = {}", total_num_rows));
    LOG_INFO(std::format("  batches = {}", batches));
}

namespace tdbs = tiledbsoma;
void test_arrow(const std::string& uri) {
    auto ctx = std::make_shared<SOMAContext>();
    const std::vector<std::string>& colnames{"n_counts", "n_genes", "louvain"};
    auto obs = tdbs::SOMAArray::open(OpenMode::read, uri, ctx);
    auto obs_mq = ManagedQuery(*obs, ctx->tiledb_ctx(), "");
    // Getting next batch:  std::optional<std::shared_ptr<ArrayBuffers>>
    auto obs_data = obs_mq.read_next();
    if (!obs_mq.results_complete()) {
        tdbs::LOG_WARN(std::format("Read of '{}' incomplete", uri));
        exit(-1);
    }
    tdbs::LOG_INFO(std::format(
        "Read complete with {} obs and {} cols",
        obs_data->get()->num_rows(),
        obs_data->get()->names().size()));
    std::vector<std::string> names = obs_data->get()->names();
    for (auto nm : names) {
        auto buf = obs_data->get()->at(nm);
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);
        ArrowSchema* schema = pp.second.get();
        tdbs::LOG_INFO(std::format(
            "Accessing '{}', retrieved '{}', n_children {}",
            nm,
            schema->name,
            schema->n_children));
    }
}

int main(int argc, char** argv) {
    LOG_CONFIG("debug");

    if (argc < 2) {
        printf("Run with CI test SOMA:\n\n");
        printf("  %s data/soco/pbmc3k_processed\n", argv[0]);
        return 0;
    }

    try {
        test_arrow(argv[1]);
        //        test_sdf(argv[1]);
    } catch (const std::exception& e) {
        printf("%s\n", e.what());
        return 1;
    }

    return 0;
};
