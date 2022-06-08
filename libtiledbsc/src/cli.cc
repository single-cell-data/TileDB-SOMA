// This file is currently a sandbox for C++ API experiments

#include <tiledbsc/tiledbsc>
#include "tiledbsc/logger_public.h"

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
    conf["config.logging_level"] = "5";
    //    Context ctx(conf);

    auto soma = SOMA::open(soma_uri, conf);
    auto array = soma.open_array("obs");
    auto ctx = soma.context();  // only need this for query, move to SOMQuery
    Query query(*ctx, array, TILEDB_READ);
    query.set_layout(TILEDB_UNORDERED);

    std::string obs_id = "obs_id";
    std::string count = "percent_mito";

    auto est = query.est_result_size_var(obs_id);
    LOG_INFO(fmt::format("est = {} {}", est[0], est[1]));

    auto buffer = ColumnBuffer::create(array, obs_id, 1 << 22);
    buffer.attach(query);
    auto buffer1 = ColumnBuffer::create(array, count, 1 << 22);
    buffer1.attach(query);

    Query::Status status;
    do {
        // Submit query and get status
        query.submit();
        status = query.query_status();
        buffer.resize(query);
        buffer1.resize(query);

        // IMPORTANT: check if there are any results, as your buffer
        // could have been too small to fit even a single result
        auto results = query.result_buffer_elements()[obs_id].first;
        bool has_results = query.result_buffer_elements()[obs_id].second != 0;
        auto offsets = buffer.offsets();

        LOG_INFO(fmt::format("status={} results={}", status, results));
        if (status == Query::Status::INCOMPLETE && !has_results) {
            // You need to reallocate your buffers, otherwise
            // you will get an infinite loop
            printf("buffers too small\n");
            return;
        } else if (has_results) {
            //            auto data = (float*)bg.data_.data();
            //            auto data = std::string_view(
            //(char*)bg.data_.data(), bg.offsets_.value()[1]);
            // auto value = "TBD";  // std::string(data.begin(), offsets[1]);
            auto value = buffer.string_view(0);
            LOG_INFO(fmt::format("value = {} len = {}", value, value.size()));
            LOG_INFO(fmt::format("value = {}", buffer1.data<float>()[0]));

            value = buffer.string_view(results - 1);
            LOG_INFO(
                fmt::format("last_value = {} len = {}", value, value.size()));
            LOG_INFO(
                fmt::format("value = {}", buffer1.data<float>()[results - 1]));
        }
    } while (status == Query::Status::INCOMPLETE);

    // Close the array
    array.close();
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
