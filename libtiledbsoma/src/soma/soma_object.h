/**
 * @file   soma_object.h
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
 * This file defines the SOMAObject class. SOMAObject is an abstract base
 * class for all public SOMA classes: SOMACollection, SOMAExperiment,
 * SOMAMeasurement, SOMA{Sparse,Dense}NdArray, and SOMADataFrame.
 */

#ifndef SOMA_OBJECT
#define SOMA_OBJECT

#include <map>
#include <string>
#include <tiledb/tiledb>

namespace tiledbsoma {

using namespace tiledb;
class SOMAObject {
   public:
    //===================================================================
    //= public non-static
    //===================================================================
    virtual ~SOMAObject() = default;

    /**
     * @brief Return a constant string describing the type of the object.
     *
     * @return std::string SOMA type
     */
    virtual std::string type() const = 0;

    /**
     * @brief Get URI of the SOMAObject.
     *
     * @return std::string URI
     */
    virtual const std::string& uri() const = 0;

    /**
     * Get the context associated with the SOMAObject.
     *
     * @return std::shared_ptr<Context>
     */
    virtual std::shared_ptr<Context> ctx() = 0;
};
}  // namespace tiledbsoma

#endif  // SOMA_OBJECT