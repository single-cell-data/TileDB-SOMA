/**
 * @file   soma_dataframe.h
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
 *   This file defines the SOMADataFrame class.
 */

#ifndef SOMA_DATAFRAME
#define SOMA_DATAFRAME

#include <map>
#include <memory>
#include <string>
#include <tiledb/tiledb>
#include "soma_object.h"

namespace tiledbsoma {

class SOMAArray;
class ArrayBuffers;

using namespace tiledb;

class SOMADataFrame : SOMAObject {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Open and return a SOMADataFrame object at the given URI.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param uri URI of the array
     * @param name Name of the array
     * @param platform_config Config parameter dictionary
     * @param column_names Columns to read
     * @param batch_size Batch size
     * @param result_order Result order
     * @param timestamp Timestamp
     * @return std::unique_ptr<SOMADataFrame> SOMADataFrame
     */
    static std::unique_ptr<SOMADataFrame> open(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::string_view name = "unnamed",
        std::map<std::string, std::string> platform_config = {},
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        std::string_view result_order = "auto",
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Open and return a SOMADataFrame object at the given URI.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param ctx TileDB context
     * @param uri URI of the array
     * @param name Name of the array
     * @param column_names Columns to read
     * @param batch_size Batch size
     * @param result_order Result order
     * @param timestamp Timestamp
     * @return std::unique_ptr<SOMADataFrame> SOMADataFrame
     */
    static std::unique_ptr<SOMAArray> open(
        tiledb_query_type_t mode,
        std::shared_ptr<Context> ctx,
        std::string_view uri,
        std::string_view name = "unnamed",
        std::vector<std::string> column_names = {},
        std::string_view batch_size = "auto",
        std::string_view result_order = "auto",
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMADataFrame object.
     *
     * @param mode TILEDB_READ or TILEDB_WRITE
     * @param uri URI of the array
     * @param name name of the array
     * @param ctx TileDB context
     * @param column_names Columns to read
     * @param batch_size Batch size
     * @param result_order Result order
     * @param timestamp Timestamp
     */
    SOMADataFrame(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::string_view name,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        std::string_view batch_size,
        std::string_view result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);
    // std::map<std::string, std::string> platform_config);

    /**
     * Closes the SOMADataFrame object.
     */
    void close();

    /**
     * Returns the constant "SOMADataFrame".
     *
     * @return std::string
     */
    std::string type() const {
        return "SOMADataFrame";
    }

    /**
     * @brief Get URI of the SOMADataFrame.
     *
     * @return std::string URI
     */
    std::string uri() const;

    /**
     * Return data schema, in the form of a TileDB ArraySchema.
     *
     * @return std::shared_ptr<ArraySchema>
     */
    std::shared_ptr<ArraySchema> schema() const;

    /**
     * Return index (dimension) column names.
     *
     * @return std::vector<std::string>
     */
    std::vector<std::string> index_column_names() const;

    /**
     * Return the number of rows in the SOMADataFrame.
     *
     * @return int64_t
     */
    int64_t count() const;

    /**
     * @brief Read the next chunk of results from the query. If all results have
     * already been read, std::nullopt is returned.
     */
    std::optional<std::shared_ptr<ArrayBuffers>> read_next();

    /**
     * @brief Write ArrayBuffers data to the dataframe.
     */
    void write(std::shared_ptr<ArrayBuffers>);

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // SOMAArray
    std::unique_ptr<SOMAArray> array_;
};
}  // namespace tiledbsoma

#endif  // SOMA_DATAFRAME