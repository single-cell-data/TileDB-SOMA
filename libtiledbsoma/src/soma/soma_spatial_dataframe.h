/**
 * @file   soma_spatial_dataframe.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2023-2024 TileDB, Inc.
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
 *   This file defines the SOMASpatialDataFrame class.
 */

#ifndef SOMA_SPATIAL_DATAFRAME
#define SOMA_SPATIAL_DATAFRAME

#include <filesystem>

#include "soma_array.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMASpatialDataFrame : public SOMAArray {
    public:
    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMADataFrame object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param column_names Columns to read
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp Timestamp
     */
    SOMASpatialDataFrame(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::vector<std::string> column_names,
        ResultOrder result_order,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMAArray(
              mode,
              uri,
              ctx,
              std::filesystem::path(uri).filename().string(),  // array name
              column_names,
              "auto",  // batch_size
              result_order,
              timestamp) {
    }

    SOMASpatialDataFrame(const SOMAArray& other)
        : SOMAArray(other) {
    }

    SOMASpatialDataFrame() = delete;
    SOMASpatialDataFrame(const SOMASpatialDataFrame&) = default;
    SOMASpatialDataFrame(SOMASpatialDataFrame&&) = delete;
    ~SOMASpatialDataFrame() = default;

    using SOMAArray::open;

    virtual void read_region() = 0;

    /**
     * Return the data schema, in the form of a ArrowSchema.
     *
     * @return std::unique_ptr<ArrowSchema>
     */
    std::unique_ptr<ArrowSchema> schema() const;

    /**
     * Return the index (dimension) column names.
     *
     * @return std::vector<std::string>
     */
    const std::vector<std::string> index_column_names() const;

    /**
     * Return the number of rows.
     *
     * @return int64_t
     */
    uint64_t count();
};
}

#endif  // SOMA_SPATIAL_DATAFRAME