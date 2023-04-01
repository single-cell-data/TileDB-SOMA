/**
 * @file   soma_group_writer.h
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
 *   This declares the SOMAGroupWriter
 */

#ifndef SOMA_GROUP_WRITER
#define SOMA_GROUP_WRITER

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <future>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>  // for group support

namespace tiledbsoma {
using namespace tiledb;

class SOMAGroupWriter {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open an array at the specified URI and return SOMAGroupWriter
     * object.
     *
     * @param uri URI of the array
     * @param platform_config Config parameter dictionary
     * @return std::unique_ptr<SOMAGroupWriter> SOMAGroupWriter
     */
    __attribute__((visibility("default"))) static void create(
        std::string_view uri,
        std::map<std::string, std::string> platform_config = {});

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMAGroupWriter object
     *
     * @param uri URI of the array
     * @param name name of the array
     * @param TODO
     */
    SOMAGroupWriter(
        std::string_view uri,
        std::map<std::string, std::string> platform_config = {},
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Construct a new SOMAGroupWriter object
     *
     * @param uri URI of the array
     * @param name name of the array
     * @param ctx TileDB context
     */
    SOMAGroupWriter(
        std::string_view uri,
        std::shared_ptr<Context> ctx,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    __attribute__((
        visibility("default"))) static std::unique_ptr<SOMAGroupWriter>
    open(
        std::string_view uri,
        std::map<std::string, std::string> platform_config = {},
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    SOMAGroupWriter() = delete;
    SOMAGroupWriter(const SOMAGroupWriter&) = delete;
    SOMAGroupWriter(SOMAGroupWriter&&) = default;
    ~SOMAGroupWriter() = default;

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // SOMAGroupWriter URI
    std::string uri_;

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // Read timestamp range (start, end)
    std::optional<std::pair<uint64_t, uint64_t>> timestamp_;
};

}  // namespace tiledbsoma

#endif  // SOMA_GROUP_WRITER
