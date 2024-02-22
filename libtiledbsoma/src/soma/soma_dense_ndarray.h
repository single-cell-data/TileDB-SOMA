/**
 * @file   soma_dense_ndarray.h
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
 *   This file defines the SOMADenseNDArray class.
 */

#ifndef SOMA_DENSE_NDARRAY
#define SOMA_DENSE_NDARRAY

#include <filesystem>

#include "soma_array.h"

namespace tiledbsoma {

class ArrayBuffers;

using namespace tiledb;

class SOMADenseNDArray : public SOMAArray {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMADenseNDArray object at the given URI.
     *
     * @param uri URI to create the SOMADenseNDArray
     * @param schema TileDB ArraySchema
     * @param platform_config Optional config parameter dictionary
     * @return std::shared_ptr<SOMADenseNDArray> opened in read mode
     */
    static std::unique_ptr<SOMADenseNDArray> create(
        std::string_view uri,
        ArraySchema schema,
        std::map<std::string, std::string> platform_config = {});

    /**
     * @brief Create a SOMADenseNDArray object at the given URI.
     *
     * @param uri URI to create the SOMADenseNDArray
     * @param schema TileDB ArraySchema
     * @param ctx TileDB context
     * @return std::shared_ptr<SOMADenseNDArray> opened in read mode
     */
    static std::unique_ptr<SOMADenseNDArray> create(
        std::string_view uri, ArraySchema schema, std::shared_ptr<Context> ctx);

    /**
     * @brief Open and return a SOMADenseNDArray object at the given URI.
     *
     * @param mode read or write
     * @param uri URI to create the SOMADenseNDArray
     * @param column_names A list of column names to use as user-defined index
     * columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must
     * exist in the schema, and at least one index column name is required.
     * @param platform_config Platform-specific options used to create this
     * SOMADenseNDArray
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @return std::shared_ptr<SOMADenseNDArray> SOMADenseNDArray
     */
    static std::unique_ptr<SOMADenseNDArray> open(
        std::string_view uri,
        OpenMode mode,
        std::map<std::string, std::string> platform_config = {},
        std::vector<std::string> column_names = {},
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Open and return a SOMADenseNDArray object at the given URI.
     *
     * @param mode read or write
     * @param ctx TileDB context
     * @param uri URI to create the SOMADenseNDArray
     * @param schema TileDB ArraySchema
     * @param column_names A list of column names to use as user-defined index
     * columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must
     * exist in the schema, and at least one index column name is required.
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp If specified, overrides the default timestamp used to
     * open this object. If unset, uses the timestamp provided by the context.
     * @return std::shared_ptr<SOMADenseNDArray> SOMADenseNDArray
     */
    static std::unique_ptr<SOMADenseNDArray> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names = {},
        ResultOrder result_order = ResultOrder::automatic,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt);

    /**
     * @brief Check if the SOMADenseNDArray exists at the URI.
     *
     * @param uri URI to create the SOMADenseNDArray
     */
    static bool exists(std::string_view uri);

    //===================================================================
    //= public non-static
    //===================================================================

    /**
     * @brief Construct a new SOMADenseNDArray object.
     *
     * @param mode read or write
     * @param uri URI of the array
     * @param ctx TileDB context
     * @param result_order Read result order: automatic (default), rowmajor, or
     * colmajor
     * @param timestamp Timestamp
     */
    SOMADenseNDArray(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<Context> ctx,
        std::vector<std::string> column_names,
        ResultOrder result_order,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp)
        : SOMAArray(
              mode,
              uri,
              std::string(std::filesystem::path(uri).filename()),  // array name
              ctx,
              column_names,
              "auto",  // batch_size
              result_order,
              timestamp) {
    }

    using SOMAArray::open;

    /**
     * Return whether the SOMADenseNDArray is sparse.
     *
     * @return false
     */
    bool is_sparse() {
        return false;
    }

    /**
     * Return the data schema, in the form of an ArrowSchema.
     *
     * @return std::unique_ptr<ArrowSchema>
     */
    std::unique_ptr<ArrowSchema> schema() const;
};
}  // namespace tiledbsoma

#endif  // SOMA_DENSE_NDARRAY
