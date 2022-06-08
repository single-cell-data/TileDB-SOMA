// This file is currently a sandbox for C++ API experiments

// TODO fixes build error on VS2019 due to "missing" include in <tiledb/type.h>
#include <stdexcept>

#include <tiledbsc/buffer_set.h>
#include <tiledbsc/logger_private.h>
#include <tiledbsc/managed_query.h>
#include <tiledbsc/query_result.h>
#include <tiledbsc/soma.h>
#include <tiledbsc/soma_collection.h>
#include <tiledbsc/tiledbsc.h>

using namespace tiledbsc;

void walk_soco(std::string_view uri) {
    auto soco = SOMACollection::open(uri);
    auto somas = soco.list_somas();

    LOG_INFO("walking soco URI = '{}'", uri);

    for (auto& [name, uri] : somas) {
        LOG_INFO("  soma {} = {}", name, uri);

        auto soma = SOMA::open(uri);
        auto arrays = soma.list_arrays();
        for (auto& [name, uri] : arrays) {
            LOG_INFO("    array {} = {}", name, uri);
        }
    }
}

void query_soma(std::string_view uri) {
    Config conf;
    conf["config.logging_level"] = "5";
    Context ctx(conf);
    auto soma = SOMA::open(uri, ctx);

    /*
        std::vector<std::string> array_names{"X/data", "obs", "var"};
        for (auto& array_name : array_names) {
            auto array = soma.open_array(array_name);
            ManagedQuery mq(array);
            auto qr = mq.execute();

            LOG_INFO("Array name = {}", array_name);
            for (auto& name : qr->names()) {
                LOG_INFO("  name = {}", name);
            }
        }
    */

    auto x_uri = soma.list_arrays()["X/data"];
    Array array(ctx, x_uri, TILEDB_READ);
    Query query(ctx, array, TILEDB_READ);
    query.set_layout(TILEDB_UNORDERED);

    std::string attr = "value";
    auto est = query.est_result_size(attr);
    printf("est = %ld\n", est);

    auto schema = array.schema();
    auto bg = BufferSet::from_attribute(schema.attribute(attr), est >> 5);

    query.set_data_buffer(
        attr, (void*)bg.data_.data(), bg.data_.size() / bg.elem_nbytes());

    Query::Status status;
    do {
        // Submit query and get status
        query.submit();
        status = query.query_status();

        // IMPORTANT: check if there are any results, as your buffer
        // could have been too small to fit even a single result
        auto results = query.result_buffer_elements()[attr].second;
        LOG_INFO("results = {}", results);
        if (status == Query::Status::INCOMPLETE && !results) {
            // You need to reallocate your buffers, otherwise
            // you will get an infinite loop
            printf("buffers too small\n");
            return;
        } else if (results) {
            auto data = (float*)bg.data_.data();
            LOG_INFO("value = {}", data[0]);
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
        LOG_FATAL("URI '{}' is not a SOMACollection. {}", argv[1], e.what());
    }

    try {
        query_soma(argv[2]);
    } catch (const std::exception& e) {
        LOG_FATAL("{}", e.what());
    }

    return 0;
};
