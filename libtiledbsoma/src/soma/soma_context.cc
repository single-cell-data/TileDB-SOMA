/**
 * @file   soma_context.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAContext class.
 */
#include <regex>
#include <thread>

#include <thread_pool/thread_pool.h>
#include "../utils/common.h"
#include "common/logging/impl/logger.h"
#include "soma_context.h"

namespace tiledbsoma {

std::shared_ptr<ThreadPool>& SOMAContext::thread_pool() {
    const std::lock_guard<std::mutex> lock(thread_pool_mutex_);
    // The first thread that gets here will create the context thread pool
    if (thread_pool_ == nullptr) {
        auto cfg = tiledb_config();
        auto concurrency = std::thread::hardware_concurrency();
        if (cfg.find(CONFIG_KEY_COMPUTE_CONCURRENCY_LEVEL) != cfg.end()) {
            auto value_str = cfg[CONFIG_KEY_COMPUTE_CONCURRENCY_LEVEL];
            try {
                concurrency = std::stoull(value_str);
            } catch (const std::exception& e) {
                throw TileDBSOMAError(
                    fmt::format(
                        "[SOMAContext] Error parsing {}: '{}' ({}) - must be a "
                        "postive integer.",
                        CONFIG_KEY_COMPUTE_CONCURRENCY_LEVEL,
                        value_str,
                        e.what()));
            }
        }

        int thread_count = std::min(std::max(1u, concurrency), 1024u);
        thread_pool_ = std::make_shared<ThreadPool>(thread_count);
    }
    return thread_pool_;
}

std::string SOMAContext::data_protocol(const std::string& uri) const {
    auto data_protocol = ctx_->data_protocol(uri);
    switch (data_protocol) {
        case tiledb::Context::DataProtocol::v2:
            return "tiledbv2";
        case tiledb::Context::DataProtocol::v3:
            return "tiledbv3";
        default:
            throw TileDBSOMAError(
                "Internal error: unrecognized TileDB data protocol. Currently only 'tiledbv2' and 'tiledbv3' are "
                "recognized.");
    }
}

void SOMAContext::validate_create_uri(const std::string_view uri) const {
    // No checks for tiledbv2. Throw error if unrecognized data protocol.
    auto data_protocol = ctx_->data_protocol(std::string(uri));
    switch (data_protocol) {
        case tiledb::Context::DataProtocol::v2:
            return;
        case tiledb::Context::DataProtocol::v3:
            break;
        default:
            throw TileDBSOMAError(
                "Internal error: unrecognized TileDB data protocol. Currently only 'tiledbv2' and 'tiledbv3' are "
                "recognized.");
    }

    std::regex storage_uri_regex("^tiledb://.*/.*://.*$", std::regex_constants::ECMAScript);
    if (std::regex_match(uri.data(), storage_uri_regex)) {
        throw TileDBSOMAError(
            fmt::format("Unsupported URI format '{}'. This format is not supported on TileDB Carrara.", uri));
    }
}

}  // namespace tiledbsoma
