/**
 * @file   soma_geometry_dataframe.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
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
#include "../utils/util.h"

#include <regex>

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMAGeometryDataFrame::create(
    std::string_view uri,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    ArrowTable spatial_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    std::vector<std::string> spatial_axes;
    auto tiledb_schema = ArrowAdapter::tiledb_schema_from_arrow_schema(
        ctx->tiledb_ctx(),
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        "SOMAGeometryDataFrame",
        true,
        platform_config,
        ArrowTable(
            std::move(spatial_columns.first),
            std::move(spatial_columns.second)));
    auto array = SOMAArray::create(
        ctx, uri, tiledb_schema, "SOMAGeometryDataFrame", timestamp);
}

std::unique_ptr<SOMAGeometryDataFrame> SOMAGeometryDataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp) {
    return std::make_unique<SOMAGeometryDataFrame>(
        mode, uri, ctx, column_names, result_order, timestamp);
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

//===================================================================
//= public non-static
//===================================================================

std::unique_ptr<ArrowSchema> SOMAGeometryDataFrame::schema() const {
    return this->arrow_schema();
}

const std::vector<std::string> SOMAGeometryDataFrame::index_column_names()
    const {
    return this->dimension_names();
}

const std::vector<std::string> SOMAGeometryDataFrame::spatial_column_names()
    const {
    std::vector<std::string> names;
    std::unordered_set<std::string> unique_names;
    std::regex rgx("tiledb__internal__(\\S+)__");
    std::smatch matches;
    for (auto dimension : this->dimension_names()) {
        if (std::regex_search(dimension, matches, rgx)) {
            if (unique_names.count(matches[1].str()) == 0) {
                unique_names.insert(matches[1].str());
                names.push_back(matches[1].str());
            }
        }
    }

    return names;
}

uint64_t SOMAGeometryDataFrame::count() {
    return this->nnz();
}

}  // namespace tiledbsoma