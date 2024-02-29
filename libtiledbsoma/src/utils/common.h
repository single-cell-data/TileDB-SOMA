/**
 * @file   common.h
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
 *   This declares the common functions in the API
 */

#ifndef TILEDBSOMA_COMMON_H
#define TILEDBSOMA_COMMON_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'
#include <string>
#include <tiledb/tiledb>

namespace tiledbsoma {

const std::string SOMA_OBJECT_TYPE_KEY = "soma_object_type";
const std::string ENCODING_VERSION_KEY = "soma_encoding_version";
const std::string ENCODING_VERSION_VAL = "1";

using MetadataValue = std::tuple<tiledb_datatype_t, uint32_t, const void*>;
enum MetadataInfo { dtype = 0, num, value };

using TimestampRange = std::pair<uint64_t, uint64_t>;

class TileDBSOMAError : public std::runtime_error {
   public:
    explicit TileDBSOMAError(const char* m)
        : std::runtime_error(m){};
    explicit TileDBSOMAError(std::string m)
        : std::runtime_error(m.c_str()){};

   public:
    virtual const char* what() const noexcept override {
        return std::runtime_error::what();
    }
};

};  // namespace tiledbsoma

#endif  // TILEDBSOMA_COMMON_H
