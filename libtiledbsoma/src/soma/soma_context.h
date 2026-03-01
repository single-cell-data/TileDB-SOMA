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
 * This file defines the SOMAContext class.
 */

#ifndef SOMA_CONTEXT
#define SOMA_CONTEXT

#include <map>
#include <memory>
#include <mutex>
#include <string>

#pragma region Forward declarations

namespace tiledb {
class Config;
class Context;
}  // namespace tiledb

namespace tiledbsoma {
class ThreadPool;
}  // namespace tiledbsoma

#pragma endregion

namespace tiledbsoma {
class SOMAContext {
    // Controls concurrency level for SOMA compute thread pool. Defaults to host
    // CPU count.
    inline static const std::string CONFIG_KEY_COMPUTE_CONCURRENCY_LEVEL = "soma.compute_concurrency_level";

   public:
    //===================================================================
    //= public non-static
    //===================================================================
    SOMAContext();

    SOMAContext(std::map<std::string, std::string> tiledb_config);

    bool operator==(const SOMAContext& other) const;

    std::shared_ptr<tiledb::Context> tiledb_ctx() const;

    std::map<std::string, std::string> tiledb_config() const;

    std::shared_ptr<ThreadPool>& thread_pool();

    /**
     * Returns the TileDB data protocol for use at a requested URI.
     *
     * @param uri The URI to get the data protocol for.
     * @returns A string description of the data protocol associated with the URI.
     */
    std::string data_protocol(const std::string& uri) const;

    /**
     * Throws a TileDBSOMA error for tiledbv3 URIs with an invalid format.
     *
     * @param uri The URI to validate.
     */
    void validate_create_uri(const std::string_view uri) const;

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<tiledb::Context> ctx_;

    // Threadpool
    std::shared_ptr<ThreadPool> thread_pool_ = nullptr;

    // Semaphore to create and use the thread_pool
    std::mutex thread_pool_mutex_;
};
}  // namespace tiledbsoma

#endif  // SOMA_CONTEXT
