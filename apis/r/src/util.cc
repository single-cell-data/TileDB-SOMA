/**
 * @file   util.cc
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
 * This file defines utilities.
 */

#include "tiledbsoma/util.h"

namespace tiledbsoma::util {

template <typename T>
VarlenBufferPair to_varlen_buffers(std::vector<T> data, bool arrow) {
    size_t nbytes = 0;
    for (auto& elem : data) {
        nbytes += elem.size();
    }

    std::vector<std::byte> result(nbytes);
    std::vector<uint64_t> offsets(data.size() + 1);
    size_t offset = 0;
    size_t idx = 0;

    for (auto& elem : data) {
        std::memcpy(result.data() + offset, elem.data(), elem.size());
        offsets[idx++] = offset;

        offset += elem.size();
    }
    offsets[idx] = offset;

    // Remove extra arrow offset when creating buffers for TileDB write
    if (!arrow) {
        offsets.pop_back();
    }

    return {result, offsets};
}

template VarlenBufferPair to_varlen_buffers(
    std::vector<std::string>, bool arrow);

bool is_tiledb_uri(std::string_view uri) {
    return uri.find("tiledb://") == 0;
}

std::string rstrip_uri(std::string_view uri) {
    return std::regex_replace(std::string(uri), std::regex("/+$"), "");
}

};  // namespace tiledbsoma::util
