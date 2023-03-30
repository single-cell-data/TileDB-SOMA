/**
 * @file   soma_group_writer.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023 TileDB, Inc.
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
 *   This file defines the SOMAGroupWriter class.
 */

#include "tiledbsoma/soma_group_writer.h"
#include "tiledbsoma/logger_public.h"
#include "tiledbsoma/util.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

// Class method
void SOMAGroupWriter::create(
    std::string_view uri, std::map<std::string, std::string> platform_config) {
    auto ctx = Context(Config(platform_config));
    tiledb::Group::create(ctx, std::string(uri));
}

std::unique_ptr<SOMAGroupWriter> SOMAGroupWriter::open(
    std::string_view uri,
    std::map<std::string, std::string> platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto ctx = std::make_shared<Context>(Config(platform_config));
    return std::make_unique<SOMAGroupWriter>(uri, ctx, timestamp);
}

//===================================================================
//= public non-static
//===================================================================

SOMAGroupWriter::SOMAGroupWriter(
    std::string_view uri,
    std::map<std::string, std::string> platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp)
    : SOMAGroupWriter(
          uri, std::make_shared<Context>(Config(platform_config)), timestamp) {
}

SOMAGroupWriter::SOMAGroupWriter(
    std::string_view uri,
    std::shared_ptr<Context> ctx,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp)
    : uri_(util::rstrip_uri(uri))
    , ctx_(ctx)
    , timestamp_(timestamp) {
    // Validate parameters
    try {
        LOG_DEBUG(fmt::format("[SOMAGroupWriter] opening group '{}'", uri_));

        // TODO: continue from here as tracked on
        // https://github.com/single-cell-data/TileDB-SOMA/issues/922

    } catch (const std::exception& e) {
        throw TileDBSOMAError(
            fmt::format("Error opening group: '{}'\n  {}", uri_, e.what()));
    }
}

}  // namespace tiledbsoma
