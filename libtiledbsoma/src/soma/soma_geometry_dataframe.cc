/**
 * @file   soma_geometry_dataframe.cc
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
 *   This file defines the SOMAGeometryDataFrame class.
 */

#include "soma_geometry_dataframe.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAGeometryDataFrame::create(
    std::string_view uri,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    std::vector<std::string> axis_names,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    
    auto tiledb_schema = ArrowAdapter::tiledb_schema_from_arrow_schema(
        ctx->tiledb_ctx(),
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        "SOMAGeometryDataFrame",
        true,
        platform_config);

    axis_names.size();
    SOMAArray::create(ctx, uri, tiledb_schema, "SOMAGeometryDataFrame", timestamp);
}

bool SOMAGeometryDataFrame::exists(
    std::string_view uri, std::shared_ptr<SOMAContext> ctx) {
    try {
        auto obj = SOMAObject::open(uri, OpenMode::read, ctx);
        return "SOMAGeometryDataFrame" == obj->type();
    } catch (TileDBSOMAError& e) {
        return false;
    }
}

void SOMAGeometryDataFrame::read_region() {
    
}
}