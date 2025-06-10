/**
 * @file   unit_soma_dense_ndarray.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMADenseNDArray class
 */

#include "common.h"

void create_soma_object(
    std::string_view soma_type,
    std::string_view uri,
    std::shared_ptr<SOMAContext> context) {
    if (soma_type == "SOMACollection") {
        return SOMACollection::create(uri, context);
    }
    if (soma_type == "SOMADataFrame") {
        auto [schema, index_columns] =
            helper::create_arrow_schema_and_index_columns(
                {helper::DimInfo(
                    {.name = "soma_joinid",
                     .tiledb_datatype = TILEDB_INT64,
                     .dim_max = 100,
                     .string_lo = "N/A",
                     .string_hi = "N/A"})},
                {helper::AttrInfo(
                    {.name = "attr1",
                     .tiledb_datatype = TILEDB_STRING_ASCII})});
        return SOMADataFrame::create(uri, schema, index_columns, context);
    }
    if (soma_type == "SOMASparseNDArray") {
        auto index_columns = helper::create_column_index_info({helper::DimInfo(
            {.name = "soma_dim_0",
             .tiledb_datatype = TILEDB_INT64,
             .dim_max = 1000,
             .string_lo = "N/A",
             .string_hi = "N/A"})});
        return SOMASparseNDArray::create(uri, "g", index_columns, context);
    }
    if (soma_type == "SOMADenseNDArray") {
        auto index_columns = helper::create_column_index_info({helper::DimInfo(
            {.name = "soma_dim_0",
             .tiledb_datatype = TILEDB_INT64,
             .dim_max = 1000,
             .string_lo = "N/A",
             .string_hi = "N/A"})});
        return SOMADenseNDArray::create(uri, "g", index_columns, context);
    }

    INFO("No support for testing type " + std::string(soma_type));
    REQUIRE(false);
}

TEST_CASE(
    "SOMAObject: basics",
    "[SOMAObject][SOMACollection][SOMADataFrame][SOMASparseNDArray]["
    "SOMADenseNDArray]") {
    auto soma_type = GENERATE(
        "SOMACollection",
        "SOMADataFrame",
        "SOMASparseNDArray",
        "SOMADenseNDArray");
    INFO("SOMA type: " + std::string(soma_type));

    std::string uri = "mem://test-soma-object-basics";
    auto context = std::make_shared<SOMAContext>();

    create_soma_object(soma_type, uri, context);

    // Close and test opening in different modes.
    OpenMode mode = GENERATE(OpenMode::read, OpenMode::write, OpenMode::del);
    INFO("Setting open mode to " + open_mode_to_string(mode));

    auto soma_obj = SOMAObject::open(uri, mode, context, std::nullopt);

    auto actual_type = soma_obj->type();
    REQUIRE(actual_type.has_value());
    REQUIRE(actual_type.value() == soma_type);

    REQUIRE(soma_obj->is_open());
    UNSCOPED_INFO(
        "Actual open mode is " + open_mode_to_string(soma_obj->mode()));
    REQUIRE(soma_obj->mode() == mode);

    soma_obj->close();
    REQUIRE(!soma_obj->is_open());
}
