/**
 * @file   soma_measurement.cc
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

#include "soma_measurement.h"
#include "soma_collection.h"
#include "soma_experiment.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMAMeasurement> SOMAMeasurement::create(
    std::shared_ptr<Context> ctx, std::string_view uri, ArraySchema schema) {
    std::string exp_uri(uri);

    SOMAGroup::create(ctx, exp_uri, "SOMAMeasurement");
    SOMADataFrame::create(ctx, exp_uri + "/var", schema);
    SOMACollection::create(ctx, exp_uri + "/X");
    SOMACollection::create(ctx, exp_uri + "/obsm");
    SOMACollection::create(ctx, exp_uri + "/obsp");
    SOMACollection::create(ctx, exp_uri + "/varm");
    SOMACollection::create(ctx, exp_uri + "/varp");

    auto group = SOMAGroup::open(TILEDB_WRITE, ctx, uri);
    group->add_member(exp_uri + "/var", false, "var");
    group->add_member(exp_uri + "/X", false, "X");
    group->add_member(exp_uri + "/obsm", false, "obsm");
    group->add_member(exp_uri + "/obsp", false, "obsp");
    group->add_member(exp_uri + "/varm", false, "varm");
    group->add_member(exp_uri + "/varp", false, "varp");
    group->close();

    return std::make_unique<SOMAMeasurement>(TILEDB_READ, uri, ctx);
}
}  // namespace tiledbsoma
