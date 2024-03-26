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
    std::string_view uri,
    ArraySchema schema,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    std::string exp_uri(uri);

    auto soma_group = SOMAGroup::create(
        ctx, exp_uri, "SOMAMeasurement", timestamp);
    SOMADataFrame::create(exp_uri + "/var", schema, ctx, timestamp);
    SOMACollection::create(exp_uri + "/X", ctx, timestamp);
    SOMACollection::create(exp_uri + "/obsm", ctx, timestamp);
    SOMACollection::create(exp_uri + "/obsp", ctx, timestamp);
    SOMACollection::create(exp_uri + "/varm", ctx, timestamp);
    SOMACollection::create(exp_uri + "/varp", ctx, timestamp);

    soma_group->set(exp_uri + "/var", URIType::absolute, "var");
    soma_group->set(exp_uri + "/X", URIType::absolute, "X");
    soma_group->set(exp_uri + "/obsm", URIType::absolute, "obsm");
    soma_group->set(exp_uri + "/obsp", URIType::absolute, "obsp");
    soma_group->set(exp_uri + "/varm", URIType::absolute, "varm");
    soma_group->set(exp_uri + "/varp", URIType::absolute, "varp");
    return std::make_unique<SOMAMeasurement>(*soma_group);
}

std::unique_ptr<SOMAMeasurement> SOMAMeasurement::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    return std::make_unique<SOMAMeasurement>(mode, uri, ctx, timestamp);
}
}  // namespace tiledbsoma
