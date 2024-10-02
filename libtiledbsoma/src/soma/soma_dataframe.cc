/**
 * @file   soma_dataframe.cc
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

#include "soma_dataframe.h"
#include "../utils/logger.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMADataFrame::create(
    std::string_view uri,
    std::unique_ptr<ArrowSchema> schema,
    ArrowTable index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    auto tiledb_schema = ArrowAdapter::tiledb_schema_from_arrow_schema(
        ctx->tiledb_ctx(),
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        "SOMADataFrame",
        true,
        platform_config);
    SOMAArray::create(ctx, uri, tiledb_schema, "SOMADataFrame", timestamp);
}

std::unique_ptr<SOMADataFrame> SOMADataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp) {
    return std::make_unique<SOMADataFrame>(
        mode, uri, ctx, column_names, result_order, timestamp);
}

std::unique_ptr<SOMADataFrame> SOMADataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::string_view name,
    std::map<std::string, std::string> platform_config) {
    return std::make_unique<SOMADataFrame>(mode, uri, name, platform_config);
}

bool SOMADataFrame::exists(
    std::string_view uri, std::shared_ptr<SOMAContext> ctx) {
    try {
        auto obj = SOMAObject::open(uri, OpenMode::read, ctx);
        return "SOMADataFrame" == obj->type();
    } catch (TileDBSOMAError& e) {
        return false;
    }
}

void SOMADataFrame::update_dataframe_schema(
    std::string uri,
    std::shared_ptr<SOMAContext> ctx,
    std::vector<std::string> drop_attrs,
    std::map<std::string, std::string> add_attrs,
    std::map<std::string, std::pair<std::string, bool>> add_enmrs) {
    ArraySchemaEvolution se(*ctx->tiledb_ctx());
    for (auto key_name : drop_attrs) {
        LOG_DEBUG(fmt::format(
            "[SOMADataFrame::update_dataframe_schema] drop col name {}",
            key_name));
        se.drop_attribute(key_name);
    }
    for (auto add_attr : add_attrs) {
        auto [attr_name, attr_type] = add_attr;

        Attribute attr(
            *ctx->tiledb_ctx(),
            attr_name,
            ArrowAdapter::to_tiledb_format(attr_type));

        if (ArrowAdapter::arrow_is_string_type(attr_type.c_str())) {
            attr.set_cell_val_num(TILEDB_VAR_NUM);
        }

        FilterList filter_list(*ctx->tiledb_ctx());
        filter_list.add_filter(Filter(*ctx->tiledb_ctx(), TILEDB_FILTER_ZSTD));
        attr.set_filter_list(filter_list);

        // An update can create (or drop) columns, or mutate existing
        // ones. A brand-new column might have nulls in it -- or it
        // might not.  And a subsequent mutator-update might set null
        // values to non-null -- or vice versa. Therefore we must be
        // careful to set nullability for all types.
        //
        // Note: this must match what DataFrame.create does:
        //
        // * DataFrame.create sets nullability for obs/var columns on
        //   initial ingest
        // * Here, we set nullability for obs/var columns on update_obs
        //
        // Users should get the same behavior either way.
        //
        // Note: this is specific to tiledbsoma.io.
        //
        // * In the SOMA API -- e.g. soma.DataFrame.create -- users
        //   bring their own Arrow schema (including nullabilities) and
        //   we must do what they say.
        // * In the tiledbsoma.io API, users bring their AnnData
        //   objects, and we compute Arrow schemas on their behalf, and
        //   we must accommodate reasonable/predictable needs.
        attr.set_nullable(true);

        // Non-enum columns:
        //
        // * add_attrs: attr_name -> Arrow type string like "i" or "U"
        // * add_enmrs: no key present
        //
        // Enum columns:
        //
        // * add_attrs: attr_name -> attr_type: Arrow type string for the index
        //   type, e.g. 'c' for int8
        // * add_enmrs: attr_name -> pair of:
        //   o enmr_type: Arrow type string the value type, e.g. "f" or "U"
        //   o ordered: bool

        auto enmr_it = add_enmrs.find(attr_name);
        bool has_enmr = enmr_it != add_enmrs.end();

        if (has_enmr) {
            auto [enmr_type, ordered] = enmr_it->second;
            LOG_DEBUG(fmt::format(
                "[SOMADataFrame::update_dataframe_schema] add col name {} "
                "index_type "
                "{} value_type {} ordered {}",
                attr_name,
                attr_type,
                enmr_type,
                ordered));
            se.add_enumeration(Enumeration::create_empty(
                *ctx->tiledb_ctx(),
                attr_name,
                ArrowAdapter::to_tiledb_format(enmr_type),
                enmr_type == "u" || enmr_type == "z" || enmr_type == "U" ||
                        enmr_type == "Z" ?
                    TILEDB_VAR_NUM :
                    1,
                ordered));
            AttributeExperimental::set_enumeration_name(
                *ctx->tiledb_ctx(), attr, attr_name);
        } else {
            LOG_DEBUG(fmt::format(
                "[SOMADataFrame::update_dataframe_schema] add col name {} type "
                "{}",
                attr_name,
                attr_type));
        }

        se.add_attribute(attr);
    }
    se.array_evolve(uri);
}

//===================================================================
//= public non-static
//===================================================================

std::unique_ptr<ArrowSchema> SOMADataFrame::schema() const {
    return this->arrow_schema();
}

const std::vector<std::string> SOMADataFrame::index_column_names() const {
    return this->dimension_names();
}

uint64_t SOMADataFrame::count() {
    return this->nnz();
}

std::optional<int64_t> SOMADataFrame::maybe_soma_joinid_shape() {
    return _maybe_soma_joinid_shape();
}

std::optional<int64_t> SOMADataFrame::maybe_soma_joinid_maxshape() {
    return _maybe_soma_joinid_maxshape();
}

}  // namespace tiledbsoma
