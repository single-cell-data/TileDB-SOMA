/**
 * @file   soma_experiment.h
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
 *   This file defines the SOMAExperiment class.
 */

#ifndef SOMA_EXPERIMENT
#define SOMA_EXPERIMENT

#include <tiledb/tiledb>

#include "soma_collection.h"
#include "soma_dataframe.h"

namespace tiledbsoma {

using namespace tiledb;
class SOMAExperiment : public SOMACollection {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAExperiment object at the given URI.
     *
     * @param uri URI to create the SOMAExperiment
     * @param schema TileDB ArraySchema
     * @param platform_config Optional config parameter dictionary
     */
    static std::unique_ptr<SOMAExperiment> create(
        std::string_view uri,
        ArraySchema schema,
        std::map<std::string, std::string> platform_config = {});

    /**
     * @brief Create a SOMAExperiment object at the given URI.
     *
     * @param uri URI to create the SOMAExperiment
     * @param schema TileDB ArraySchema
     * @param ctx TileDB context
     */
    static std::unique_ptr<SOMAExperiment> create(
        std::string_view uri, ArraySchema schema, std::shared_ptr<Context> ctx);

    //===================================================================
    //= public non-static
    //===================================================================

    SOMAExperiment(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<Context> ctx,
        std::optional<std::pair<uint64_t, uint64_t>> timestamp = std::nullopt)
        : SOMACollection(mode, uri, ctx, timestamp) {
    }

    SOMAExperiment() = delete;
    SOMAExperiment(const SOMAExperiment&) = delete;
    SOMAExperiment(SOMAExperiment&&) = default;
    ~SOMAExperiment() = default;

    /**
     * Return the constant "SOMAExperiment".
     *
     * @return std::string
     */
    const std::string type() const {
        return "SOMAExperiment";
    }

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // Primary annotations on the observation axis
    std::shared_ptr<SOMADataFrame> obs_;

    // A collection of named measurements
    std::shared_ptr<SOMACollection> ms_;
};
}  // namespace tiledbsoma

#endif  // SOMA_EXPERIMENT