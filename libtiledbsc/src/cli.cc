// This file is currently a sandbox for C++ API experiments

#include <tiledbsc/tiledbsc>

using namespace tiledbsc;

void walk_soco(std::string_view uri) {
    auto soco = SOMACollection::open(uri);
    auto somas = soco.list_somas();

    LOG_INFO(fmt::format("walking soco URI = {}", uri));

    for (auto& [name, uri] : somas) {
        LOG_INFO(fmt::format("  soma {} = {}", name, uri));

        auto soma = SOMA::open(uri);
        auto arrays = soma.list_arrays();
        for (auto& [name, uri] : arrays) {
            LOG_INFO(fmt::format("    array {} = {}", name, uri));
        }
    }
}

void slice_soma(std::string_view soma_uri) {
    Config conf;
    // conf["config.logging_level"] = "5";

    auto soma = SOMA::open(soma_uri, conf);
    auto mq = ManagedQuery(soma.open_array("obs"));

    mq.select_columns({"obs_id", "percent_mito"});
    mq.select_points<std::string>("obs_id", {"AAACATACAACCAC-1"});
    mq.select_ranges<std::string>(
        "obs_id", {{"TTTCGAACTCTCAT-1", "TTTGCATGCCTCAC-1"}});

    while (auto num_cells = mq.execute()) {
        for (size_t i = 0; i < num_cells; i++) {
            LOG_INFO(fmt::format(
                "obs_id = {} percent_mito = {}",
                mq.string_view("obs_id", i),
                mq.data<float>("percent_mito")[i]));
        }
    }
}

int main(int argc, char** argv) {
    if (argc != 3) {
        printf("Usage: %s soco_uri soma_uri\n", argv[0]);
        return 1;
    }

    LOG_CONFIG("debug");

    try {
        walk_soco(argv[1]);
    } catch (const std::exception& e) {
        LOG_FATAL(fmt::format(
            "URI '{}' is not a SOMACollection. {}", argv[1], e.what()));
    }

    try {
        slice_soma(argv[2]);
    } catch (const std::exception& e) {
        LOG_FATAL(e.what());
    }

    return 0;
};
