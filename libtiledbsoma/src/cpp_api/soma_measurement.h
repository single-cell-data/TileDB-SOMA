/**
 * @file   soma_measurement.h
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
 *   This file defines the SOMAMeasurement class.
 */

#ifndef SOMA_MEASUREMENT
#define SOMA_MEASUREMENT

#include <map>
#include <memory>
#include <string>
#include <tiledb/tiledb>
#include "soma_collection.h"
#include "soma_dataframe.h"
#include "soma_dense_ndarray.h"
#include "soma_object.h"
#include "soma_sparse_ndarray.h"

namespace tiledbsoma {

class SOMAGroup;

using namespace tiledb;

class SOMAMeasurement : SOMACollection {
   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // TileDB context
    std::shared_ptr<Context> ctx_;

    // SOMAGroup
    std::unique_ptr<SOMAGroup> group_;

    // Primary annotations on the variable axis
    SOMADataFrame var_;

    // A collection of matrices, each containing measured feature values
    SOMACollection X_;

    // A collection of dense matrices containing annotations of each obs row
    SOMACollection obsm_;

    // A collection of sparse matrices containing pairwise annotations of each
    // obs row
    SOMACollection obsp_;

    // A collection of dense matrices containing annotations of each var row
    SOMACollection varm_;

    // A collection of sparse matrices containing pairwise annotations of each
    // var row
    SOMACollection varp_;
};
}  // namespace tiledbsoma

#endif  // SOMA_MEASUREMENT