// TODO fixes build error on VS2019 due to "missing" include in <tiledb/type.h>
#include <stdexcept>

#include <tiledbsc/logger_private.h>
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

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Usage: %s uri\n", argv[0]);
        return 1;
    }

    LOG_CONFIG("debug");

    try {
        walk_soco(argv[1]);
    } catch (const std::exception& e) {
        LOG_FATAL("URI '{}' is not a SOMACollection. {}", argv[1], e.what());
    }

    return 0;
};
