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
    //= public non-static
    //===================================================================

    SOMAExperiment(
        tiledb_query_type_t mode,
        std::string_view uri,
        std::shared_ptr<Context> ctx,
        SOMADataFrame& obs,
        SOMACollection& ms);

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // Primary annotations on the observation axis
    SOMADataFrame obs_;

    // A collection of named measurements
    SOMACollection ms_;
};
}  // namespace tiledbsoma

#endif  // SOMA_EXPERIMENT