/**
 * @file   util.h
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
 * This file defines the utility functions
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <regex>
#include <tiledbsoma/tiledbsoma>

namespace tiledbsoma::util {

using VarlenBufferPair =
    std::pair<std::vector<std::byte>, std::vector<uint64_t>>;

template <typename T>
VarlenBufferPair to_varlen_buffers(std::vector<T> data, bool arrow = true);

template <class T>
std::vector<T> to_vector(const tcb::span<T>& s) {
    return std::vector<T>(s.begin(), s.end());
}

/**
 * @brief Check if the provided URI is a TileDB Cloud URI.
 *
 * @param uri URI to check
 * @return true URI is a TileBD Cloud URI
 * @return false URI is not a TileBD Cloud URI
 */
bool is_tiledb_uri(std::string_view uri);

/**
 * @brief Remove all trailing '/' from URI.
 *
 * @param uri URI
 * @return std::string URI without trailing '/'
 */
std::string rstrip_uri(std::string_view uri);

}  // namespace tiledbsoma::util

#endif
