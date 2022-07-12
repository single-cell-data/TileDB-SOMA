// This file is currently a sandbox for C++ API experiments

#include <tiledbsc/tiledbsc>

using namespace tiledbsc;

void walk_soco(std::string_view uri) {
    auto soco = SOMACollection::open(uri);
    auto somas = soco->list_somas();

    LOG_INFO(fmt::format("walking soco URI = {}", uri));

    for (auto& [name, uri] : somas) {
        LOG_INFO(fmt::format("  soma {} = {}", name, uri));

        auto soma = SOMA::open(uri);
        auto arrays = soma->list_arrays();
        for (auto& [name, uri] : arrays) {
            LOG_INFO(fmt::format("    array {} = {}", name, uri));
        }
    }
}

void slice_soma(std::string_view soma_uri) {
    Config conf;
    conf["config.logging_level"] = "5";

    auto soma = SOMA::open(soma_uri, conf);
    auto array = soma->open_array("obs");
    auto mq = ManagedQuery(array);

    auto q = std::make_unique<Query>(array->schema().context(), *array);
    auto est_bytes = q->est_result_size_var("obs_id");
    LOG_DEBUG(fmt::format("est_num_cells = {}", est_bytes[0] / 8));
    LOG_DEBUG(fmt::format("est_bytes = {}", est_bytes[1]));

    mq.select_columns({"obs_id", "percent_mito"});
    mq.select_points<std::string>("obs_id", {"AAACATACAACCAC-1"});
    mq.select_ranges<std::string>(
        "obs_id", {{"TTTCGAACTCTCAT-1", "TTTGCATGCCTCAC-1"}});

    size_t total_cells = 0;
    while (!mq.is_complete()) {
        auto num_cells = mq.submit();
        LOG_DEBUG(fmt::format("num_cells = {}", num_cells));
        total_cells += num_cells;
        continue;

        auto mito = mq.data<float>("percent_mito");
        for (size_t i = 0; i < num_cells; i++) {
            auto obs = mq.string_view("obs_id", i);
            LOG_INFO(
                fmt::format("obs_id = {} percent_mito = {}", obs, mito[i]));
        }
    }
    LOG_DEBUG(fmt::format("total_cells = {}", total_cells));
}

void soma_query(std::string_view soma_uri) {
    Config conf;
    conf["config.logging_level"] = "5";

    auto soma = SOMA::open(soma_uri, conf);
    auto sq = soma->query();
    auto ctx = soma->context();

    std::string obs_attr = "cell_type";
    std::string obs_val = "B cell";
    auto obs_qc = QueryCondition::create(*ctx, obs_attr, obs_val, TILEDB_EQ);
    std::vector<std::string> obs_cols = {obs_attr};
    sq->set_obs_condition(obs_qc);
    sq->select_obs_attrs(obs_cols);

    uint64_t var_val = 50;
    auto var_qc = QueryCondition::create<uint64_t>(
        *ctx, "n_cells", var_val, TILEDB_LT);
    std::vector<std::string> var_cols = {"var_id"};
    //    sq->set_var_condition(var_qc);
    sq->select_var_attrs(var_cols);

    std::vector<std::string> obs_ids = {
        "AAACATACAACCAC-1", "AAACATTGATCAGC-1", "TTTGCATGCCTCAC-1"};
    std::vector<std::string> var_ids = {"AAGAB", "AAR2", "ZRANB3"};
    // sq->select_obs_ids(obs_ids);
    // sq->select_var_ids(var_ids);

    size_t total_cells = 0;
    while (auto results = sq->next_results()) {
        auto num_cells = results->at("obs_id")->size();
        total_cells += num_cells;
        LOG_DEBUG(fmt::format("num_cells = {}", num_cells));
        if (num_cells < 20) {
            for (size_t i = 0; i < num_cells; i++) {
                LOG_DEBUG(fmt::format(
                    "{} {} {}",
                    results->at("obs_id")->string_view(i),
                    results->at("var_id")->string_view(i),
                    results->at("value")->data<float>()[i]));
            }
        }
    }
    LOG_DEBUG(fmt::format("total_cells = {}", total_cells));
}

void soco_query(std::string_view soco_uri) {
    Config conf;
    conf["config.logging_level"] = "5";

    auto soco = SOMACollection::open(soco_uri, conf);
    SOMACollectionQuery sqs(soco.get());

    LOG_DEBUG("Submit");
    sqs.next_results();
    LOG_DEBUG("Submit again");
    sqs.next_results();
    LOG_DEBUG("Submit again");
    sqs.next_results();
    LOG_DEBUG("Done");
}

int main(int argc, char** argv) {
    if (argc != 3) {
        printf("Usage: %s soco_uri soma_uri\n", argv[0]);
        return 1;
    }

    LOG_CONFIG("debug");

    try {
        walk_soco(argv[1]);
        // slice_soma(argv[2]);
        soco_query(argv[1]);
    } catch (const std::exception& e) {
        LOG_FATAL(e.what());
    }

    return 0;
};
