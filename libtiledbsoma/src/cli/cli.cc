/**
 * @file   cli.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file is currently a sandbox for C++ API experiments
 */

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
    auto obs = SOMAArray::open(OpenMode::read, uri + "/obs", ctx, "obs");
    auto obs_data = obs->read_next();

    // Read all values from the var array
    auto var = SOMAArray::open(OpenMode::read, uri + "/ms/RNA/var", ctx, "var");
    auto var_data = var->read_next();

    // Check if obs and var reads are complete
    if (obs->results_complete() && var->results_complete()) {
        LOG_INFO("var and obs queries are complete");
    }

    // Read all values from the X/data array
    auto x_data = SOMAArray::open(
        OpenMode::read,
        uri + "/ms/RNA/X/data",
        std::make_shared<SOMAContext>(config),
        "X/data");

    int batches = 0;
    int total_num_rows = 0;

    // Handle incomplete queries
    while (auto batch = x_data->read_next()) {
        batches++;
        total_num_rows += (*batch)->num_rows();
    }

    LOG_INFO(fmt::format("X/data rows = {}", total_num_rows));
    LOG_INFO(fmt::format("  batches = {}", batches));
}

namespace tdbs = tiledbsoma;
void test_arrow(const std::string& uri) {
    const std::vector<std::string>& colnames{"n_counts", "n_genes", "louvain"};
    auto obs = tdbs::SOMAArray::open(
        OpenMode::read, uri, std::make_shared<SOMAContext>(), "", colnames);
    // Getting next batch:  std::optional<std::shared_ptr<ArrayBuffers>>
    auto obs_data = obs->read_next();
    if (!obs->results_complete()) {
        tdbs::LOG_WARN(fmt::format("Read of '{}' incomplete", uri));
        exit(-1);
    }
    tdbs::LOG_INFO(fmt::format(
        "Read complete with {} obs and {} cols",
        obs_data->get()->num_rows(),
        obs_data->get()->names().size()));
    std::vector<std::string> names = obs_data->get()->names();
    for (auto nm : names) {
        auto buf = obs_data->get()->at(nm);
        auto pp = tdbs::ArrowAdapter::to_arrow(buf);
        ArrowSchema* schema = pp.second.get();
        tdbs::LOG_INFO(fmt::format(
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
        printf("  %s test/soco/pbmc3k_processed\n", argv[0]);
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
