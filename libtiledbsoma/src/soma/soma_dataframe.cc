/**
 * @file   soma_dataframe.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMADataFrame class.
 */

#include "soma_dataframe.h"
#include <format>
#include "../utils/logger.h"
#include "soma_coordinates.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMADataFrame::create(
    std::string_view uri,
    const std::unique_ptr<ArrowSchema>& schema,
    const ArrowTable& index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    auto [tiledb_schema, soma_schema_extension] =
        ArrowAdapter::tiledb_schema_from_arrow_schema(
            ctx->tiledb_ctx(),
            schema,
            index_columns,
            std::nullopt,
            "SOMADataFrame",
            true,
            platform_config);
    SOMAArray::create(
        ctx, uri, tiledb_schema, "SOMADataFrame", std::nullopt, timestamp);
}

std::unique_ptr<SOMADataFrame> SOMADataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp) {
    auto array = std::make_unique<SOMADataFrame>(mode, uri, ctx, timestamp);

    if (!array->check_type("SOMADataFrame")) {
        throw TileDBSOMAError(
            "[SOMADataFrame::open] Object is not a SOMADataFrame");
    }

    return array;
}

std::unique_ptr<SOMADataFrame> SOMADataFrame::open(
    std::string_view uri,
    OpenMode mode,
    std::map<std::string, std::string> platform_config,
    std::optional<TimestampRange> timestamp) {
    auto array = std::make_unique<SOMADataFrame>(
        mode, uri, platform_config, timestamp);

    if (!array->check_type("SOMADataFrame")) {
        throw TileDBSOMAError(
            "[SOMADataFrame::open] Object is not a SOMADataFrame");
    }

    return array;
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
    std::vector<std::string> drop_attrs,
    std::map<std::string, std::string> add_attrs,
    std::map<std::string, std::pair<std::string, bool>> add_enmrs) {
    const auto& tctx = *ctx_->tiledb_ctx();
    ArraySchemaEvolution se(tctx);
    for (auto key_name : drop_attrs) {
        LOG_DEBUG(std::format(
            "[SOMADataFrame::update_dataframe_schema] drop col name {}",
            key_name));
        se.drop_attribute(key_name);
    }
    for (auto add_attr : add_attrs) {
        auto [attr_name, attr_type] = add_attr;

        Attribute attr(
            tctx, attr_name, ArrowAdapter::to_tiledb_format(attr_type));

        if (ArrowAdapter::arrow_is_var_length_type(attr_type.c_str())) {
            attr.set_cell_val_num(TILEDB_VAR_NUM);
        }

        FilterList filter_list(tctx);
        filter_list.add_filter(Filter(tctx, TILEDB_FILTER_ZSTD));
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
        bool request_has_enmr = enmr_it != add_enmrs.end();
        bool array_already_has_it = false;

        if (request_has_enmr) {
            auto [enmr_type, ordered] = enmr_it->second;

            // SOMA converts string and binary to large string and large binary
            auto enmr_name = attr_name + "_" +
                             ((enmr_type == "u") ?
                                  "U" :
                                  (enmr_type == "z" ? "Z" : enmr_type));

            LOG_DEBUG(std::format(
                "[SOMADataFrame::update_dataframe_schema] add col name {} "
                "index_type "
                "{} value_type {} ordered {}",
                attr_name,
                attr_type,
                enmr_type,
                ordered));

            // First we need to check if the array has (ever had) an enumeration
            // of this type in its history.
            //
            // Note that in core, attrs have names and enumerations have
            // names, and the names need not match -- and multiple attrs can
            // use the same enumeration. For tiledbsoma, though, we don't take
            // advantage of this flexibility: we use the same name for an attr
            // and its enumeration (if it has one).
            //
            // It's quite possible the array had this name as an enumerated
            // column, and then that column was dropped, and it's now being
            // re-added.
            //
            // * We need to check if this is the case.
            // * If it is, we need to check if the enumeration is compatible.
            try {
                const auto existing_enmr = ArrayExperimental::get_enumeration(
                    tctx, *arr_, enmr_name);
                auto existing_type = ArrowAdapter::tdb_to_arrow_type(
                    existing_enmr.type());
                auto existing_ordered = existing_enmr.ordered();

                // It's there but it'll cause a problem further down in core.
                // We cannot continue without intervention, but let's at least
                // be as clear as possible why we can't continue.
                if (ordered != existing_ordered || enmr_type != existing_type) {
                    throw TileDBSOMAError(std::format(
                        "[SOMADataFrame::update_dataframe_schema] requested "
                        "enumeration [type='{}', ordered={}] incompatible with "
                        "existing [type='{}', ordered={}]",
                        enmr_type,
                        ordered,
                        existing_type,
                        existing_ordered));
                }

                array_already_has_it = true;
            } catch (const std::exception& e) {
                // Oddly, catch (tiledb::TileDBError& e) did not enter this
                // block in unit test.
                //
                // Make sure the exception was for the reason we're considering.
                // If it was for some other reason, fail.
                const std::string msg = e.what();
                if (msg.find("already exists in this ArraySchema") !=
                    std::string::npos) {
                    array_already_has_it = true;
                } else if (
                    msg.find("Unable to check if unknown enumeration is "
                             "loaded. No enumeration named") !=
                    std::string::npos) {
                    array_already_has_it = false;
                } else {
                    throw(e);
                }
            }

            if (!array_already_has_it) {
                se.add_enumeration(Enumeration::create_empty(
                    tctx,
                    enmr_name,
                    ArrowAdapter::to_tiledb_format(enmr_type),
                    enmr_type == "u" || enmr_type == "z" || enmr_type == "U" ||
                            enmr_type == "Z" ?
                        TILEDB_VAR_NUM :
                        1,
                    ordered));
                AttributeExperimental::set_enumeration_name(
                    tctx, attr, enmr_name);
            }

        } else {
            LOG_DEBUG(std::format(
                "[SOMADataFrame::update_dataframe_schema] add col name {} type "
                "{}",
                attr_name,
                attr_type));
        }

        se.add_attribute(attr);
    }

    se.array_evolve(uri_);

    // When we evolve the schema, the ArraySchema needs to be updated to the
    // latest version so re-open the Array. tiledb::Array::open will reopen
    // the array at the set timestamp
    auto mode = arr_->query_type();
    arr_->close();
    arr_->open(mode);
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
