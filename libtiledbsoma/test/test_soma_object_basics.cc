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

void create_soma_object(std::string_view soma_type, std::string_view uri, std::shared_ptr<SOMAContext> context) {
    if (soma_type == "SOMACollection") {
        return SOMACollection::create(uri, context);
    }
    if (soma_type == "SOMAExperiment") {
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(
            {helper::DimInfo(
                {.name = "soma_joinid",
                 .tiledb_datatype = TILEDB_INT64,
                 .dim_max = 100,
                 .string_lo = "N/A",
                 .string_hi = "N/A"})},
            {helper::AttrInfo({.name = "attr1", .tiledb_datatype = TILEDB_STRING_ASCII})});
        return SOMAExperiment::create(uri, schema, index_columns, context);
    }
    if (soma_type == "SOMAMeasurement") {
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(
            {helper::DimInfo(
                {.name = "soma_joinid",
                 .tiledb_datatype = TILEDB_INT64,
                 .dim_max = 100,
                 .string_lo = "N/A",
                 .string_hi = "N/A"})},
            {helper::AttrInfo({.name = "attr1", .tiledb_datatype = TILEDB_STRING_ASCII})});
        return SOMAMeasurement::create(uri, schema, index_columns, context);
    }
    if (soma_type == "SOMAScene") {
        return SOMAScene::create(uri, context, std::nullopt);
    }
    if (soma_type == "SOMADataFrame") {
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(
            {helper::DimInfo(
                {.name = "soma_joinid",
                 .tiledb_datatype = TILEDB_INT64,
                 .dim_max = 100,
                 .string_lo = "N/A",
                 .string_hi = "N/A"})},
            {helper::AttrInfo({.name = "attr1", .tiledb_datatype = TILEDB_STRING_ASCII})});
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
    if (soma_type == "SOMAMultiscaleImage") {
        SOMACoordinateSpace coord_space{};
        return SOMAMultiscaleImage::create(uri, context, coord_space, std::nullopt);
    }
    if (soma_type == "SOMAPointCloudDataFrame") {
        std::vector<helper::DimInfo> dim_infos({
            helper::DimInfo(
                {.name = "soma_joinid",
                 .tiledb_datatype = TILEDB_INT64,
                 .dim_max = 1000,
                 .string_lo = "N/A",
                 .string_hi = "N/A"}),
            helper::DimInfo(
                {.name = "x",
                 .tiledb_datatype = TILEDB_UINT32,
                 .dim_max = 100,
                 .string_lo = "N/A",
                 .string_hi = "N/A"}),
            helper::DimInfo(
                {.name = "y",
                 .tiledb_datatype = TILEDB_UINT32,
                 .dim_max = 100,
                 .string_lo = "N/A",
                 .string_hi = "N/A"}),
        });
        std::vector<helper::AttrInfo> attr_infos(
            {helper::AttrInfo({.name = "radius", .tiledb_datatype = TILEDB_FLOAT64})});
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(dim_infos, attr_infos);
        SOMACoordinateSpace coord_space{};
        return SOMAPointCloudDataFrame::create(uri, schema, index_columns, coord_space, context);
    }
    if (soma_type == "SOMAGeometryDataFrame") {
        std::vector<helper::DimInfo> dim_infos(
            {helper::DimInfo(
                 {.name = "soma_joinid",
                  .tiledb_datatype = TILEDB_INT64,
                  .dim_max = 1000,
                  .string_lo = "N/A",
                  .string_hi = "N/A"}),
             helper::DimInfo(
                 {.name = "soma_geometry",
                  .tiledb_datatype = TILEDB_GEOM_WKB,
                  .dim_max = 100,
                  .string_lo = "N/A",
                  .string_hi = "N/A"})});

        std::vector<helper::AttrInfo> attr_infos(
            {helper::AttrInfo({.name = "quality", .tiledb_datatype = TILEDB_FLOAT64})});
        SOMACoordinateSpace coord_space{};
        auto [schema, index_columns] = helper::create_arrow_schema_and_index_columns(
            dim_infos, attr_infos, coord_space);
        return SOMAGeometryDataFrame::create(uri, schema, index_columns, coord_space, context);
    }

    INFO("No support for testing type " + std::string(soma_type));
    REQUIRE(false);
}

TEST_CASE(
    "SOMAObject: basics",
    "[SOMAObject][SOMACollection][SOMAMeasurement][SOMAExperiment][SOMADataFrame][SOMASparseNDArray][SOMADenseNDArray]["
    "SOMAScene][SOMAPointCloudDataFrame][SOMAMultiscaleImage]") {
    auto soma_type = GENERATE(
        "SOMACollection",
        "SOMAMeasurement",
        "SOMAExperiment",
        "SOMADataFrame",
        "SOMASparseNDArray",
        "SOMADenseNDArray",
        "SOMAScene",
        "SOMAPointCloudDataFrame",
        "SOMAGeometryDataFrame",
        "SOMAMultiscaleImage");
    INFO("SOMA type: " + std::string(soma_type));

    std::string uri = "mem://test-soma-object-basics";
    auto context = std::make_shared<SOMAContext>();

    create_soma_object(soma_type, uri, context);

    // Close and test opening in different modes.
    OpenMode mode = GENERATE(OpenMode::soma_read, OpenMode::soma_write, OpenMode::soma_delete);
    INFO("Setting open mode to " + open_mode_to_string(mode));

    auto soma_obj = SOMAObject::open(uri, mode, context, std::nullopt);

    auto actual_type = soma_obj->type();
    REQUIRE(actual_type.has_value());
    REQUIRE(actual_type.value() == soma_type);

    auto actual_type2 = get_soma_type_metadata_value(uri, *context);
    REQUIRE(actual_type == soma_type);

    REQUIRE(soma_obj->is_open());
    UNSCOPED_INFO("Actual open mode is " + open_mode_to_string(soma_obj->mode()));
    REQUIRE(soma_obj->mode() == mode);

    soma_obj->close();
    REQUIRE(!soma_obj->is_open());
}

TEST_CASE("SOMAArray: Invalid datatype metadata", "[SOMAArray][Metadata]") {
    std::string uri = "mem://test-get-array-type";

    // Create TileDB array with no metadata.
    SOMAContext ctx{};
    const auto& tiledb_ctx = *ctx.tiledb_ctx();
    Domain domain(tiledb_ctx);
    domain.add_dimension(Dimension::create<int64_t>(tiledb_ctx, "dim1", {{1, 4}}, 4));
    ArraySchema schema{tiledb_ctx, TILEDB_SPARSE};
    schema.set_domain(domain);
    schema.add_attribute(Attribute::create<int64_t>(tiledb_ctx, "attr1"));
    Array::create(tiledb_ctx, uri, schema);

    // Check SOMA type calls fail when no type metadata.
    REQUIRE_THROWS(
        get_soma_type_metadata_value(uri, ctx),
        Catch::Matchers::ContainsSubstring("Missing the required metadata key"));
    REQUIRE_THROWS(
        get_soma_type_metadata_value_from_array(uri, ctx),
        Catch::Matchers::ContainsSubstring("Missing the required metadata key"));
    REQUIRE_THROWS(
        get_soma_type_metadata_value_from_group(uri, ctx), Catch::Matchers::EndsWith("Group does not exist."));

    SECTION("Valid metadata string that is not a SOMA type") {
        // Add soma type metdata with invalid type.
        Array array{tiledb_ctx, uri, TILEDB_WRITE, TemporalPolicy()};
        std::string fake_type{"not_a_soma_type"};
        array.put_metadata(
            "soma_object_type", TILEDB_STRING_UTF8, static_cast<uint32_t>(fake_type.length()), fake_type.data());
        array.close();

        auto actual1 = get_soma_type_metadata_value(uri, ctx);
        REQUIRE(actual1 == fake_type);

        auto actual2 = get_soma_type_metadata_value_from_array(uri, ctx);
        REQUIRE(actual1 == fake_type);

        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_group(uri, ctx), Catch::Matchers::EndsWith("Group does not exist."));
    }
    SECTION("Invalid metadata datatype: wrong string type") {
        // Add soma type metdata with invalid type.
        Array array{tiledb_ctx, uri, TILEDB_WRITE, TemporalPolicy()};
        std::string fake_type{"not_a_soma_type"};
        array.put_metadata(
            "soma_object_type", TILEDB_STRING_ASCII, static_cast<uint32_t>(fake_type.length()), fake_type.data());
        array.close();

        REQUIRE_THROWS(
            get_soma_type_metadata_value(uri, ctx),
            Catch::Matchers::EndsWith("has datatype 'STRING_ASCII'. Expected datatype 'STRING_UTF8'."));
        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_array(uri, ctx),
            Catch::Matchers::EndsWith("has datatype 'STRING_ASCII'. Expected datatype 'STRING_UTF8'."));
        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_group(uri, ctx), Catch::Matchers::EndsWith("Group does not exist."));
    }
    SECTION("Invalid metadata datatype: numeric type") {
        // Add soma type metdata with invalid type.
        Array array{tiledb_ctx, uri, TILEDB_WRITE, TemporalPolicy()};
        std::string fake_type{"not_a_soma_type"};
        int64_t value{4};
        array.put_metadata("soma_object_type", TILEDB_INT64, 1, &value);
        array.close();

        REQUIRE_THROWS(
            get_soma_type_metadata_value(uri, ctx),
            Catch::Matchers::EndsWith("has datatype 'INT64'. Expected datatype 'STRING_UTF8'."));
        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_array(uri, ctx),
            Catch::Matchers::EndsWith("has datatype 'INT64'. Expected datatype 'STRING_UTF8'."));
        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_group(uri, ctx), Catch::Matchers::EndsWith("Group does not exist."));
    }
}

TEST_CASE("SOMAGroup: Invalid datatype metadata", "[SOMAGroup][Metadata]") {
    std::string uri = "mem://test-get-group-type";

    // Create TileDB group with no metadata;
    SOMAContext ctx{};
    const auto& tiledb_ctx = *ctx.tiledb_ctx();
    Group::create(tiledb_ctx, uri);

    // Check SOMA type calls fail when no type metadata.
    REQUIRE_THROWS(
        get_soma_type_metadata_value(uri, ctx),
        Catch::Matchers::ContainsSubstring("Missing the required metadata key"));
    REQUIRE_THROWS(
        get_soma_type_metadata_value_from_group(uri, ctx),
        Catch::Matchers::ContainsSubstring("Missing the required metadata key"));
    REQUIRE_THROWS(
        get_soma_type_metadata_value_from_array(uri, ctx), Catch::Matchers::EndsWith("ZArray does not exist."));

    SECTION("Valid metadata string that is not a SOMA type") {
        // Add soma type metdata with invalid type.
        Group group{tiledb_ctx, uri, TILEDB_WRITE};
        std::string fake_type{"not_a_soma_type"};
        group.put_metadata(
            "soma_object_type", TILEDB_STRING_UTF8, static_cast<uint32_t>(fake_type.length()), fake_type.data());
        group.close();

        auto actual1 = get_soma_type_metadata_value(uri, ctx);
        REQUIRE(actual1 == fake_type);

        auto actual2 = get_soma_type_metadata_value_from_group(uri, ctx);
        REQUIRE(actual1 == fake_type);

        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_array(uri, ctx), Catch::Matchers::EndsWith("Array does not exist."));
    }
    SECTION("Invalid metadata datatype: wrong string type") {
        // Add soma type metdata with invalid type.
        Group group{tiledb_ctx, uri, TILEDB_WRITE};
        std::string fake_type{"not_a_soma_type"};
        group.put_metadata(
            "soma_object_type", TILEDB_STRING_ASCII, static_cast<uint32_t>(fake_type.length()), fake_type.data());
        group.close();

        REQUIRE_THROWS(
            get_soma_type_metadata_value(uri, ctx),
            Catch::Matchers::EndsWith("has datatype 'STRING_ASCII'. Expected datatype 'STRING_UTF8'."));
        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_group(uri, ctx),
            Catch::Matchers::EndsWith("has datatype 'STRING_ASCII'. Expected datatype 'STRING_UTF8'."));
        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_array(uri, ctx), Catch::Matchers::EndsWith("Array does not exist."));
    }
    SECTION("Invalid metadata datatype: numeric type") {
        // Add soma type metdata with invalid type.
        Group group{tiledb_ctx, uri, TILEDB_WRITE};
        std::string fake_type{"not_a_soma_type"};
        int64_t value{4};
        group.put_metadata("soma_object_type", TILEDB_INT64, 1, &value);
        group.close();

        REQUIRE_THROWS(
            get_soma_type_metadata_value(uri, ctx),
            Catch::Matchers::EndsWith("has datatype 'INT64'. Expected datatype 'STRING_UTF8'."));
        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_group(uri, ctx),
            Catch::Matchers::EndsWith("has datatype 'INT64'. Expected datatype 'STRING_UTF8'."));
        REQUIRE_THROWS(
            get_soma_type_metadata_value_from_array(uri, ctx), Catch::Matchers::EndsWith("Array does not exist."));
    }
}
