/** * @file   unit_soma_geometry_dataframe.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMAGeometryDataFrame class
 */

#include <format>
#include <vector>
#include "../src/geometry/geometry.h"
#include "../src/geometry/operators/io/write.h"
#include "../src/utils/common.h"
#include "common.h"

const int64_t SOMA_JOINID_DIM_MAX = 99;
const SOMACoordinateSpace coord_space(
    {SOMAAxis{"x", std::nullopt}, SOMAAxis{"y", std::nullopt}});

TEST_CASE("SOMAGeometryDataFrame: basic", "[SOMAGeometryDataFrame]") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri{"mem://unit-test-geometry-basic"};
    PlatformConfig platform_config{};

    std::vector<helper::DimInfo> dim_infos(
        {helper::DimInfo(
             {.name = "soma_joinid",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = SOMA_JOINID_DIM_MAX,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "soma_geometry",
              .tiledb_datatype = TILEDB_GEOM_WKB,
              .dim_max = 100,
              .string_lo = "N/A",
              .string_hi = "N/A"})});

    std::vector<helper::AttrInfo> attr_infos({helper::AttrInfo(
        {.name = "quality", .tiledb_datatype = TILEDB_FLOAT64})});

    // Check the geometry dataframe doesn't exist yet.
    REQUIRE(!SOMAGeometryDataFrame::exists(uri, ctx));

    // Create the geometry dataframe.
    auto [schema, index_columns] =
        helper::create_arrow_schema_and_index_columns(
            dim_infos, attr_infos, coord_space);

    SOMAGeometryDataFrame::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        coord_space,
        ctx,
        platform_config,
        std::nullopt);

    // Check the geometry dataframe exists and it cannot be read as a
    // different object.
    REQUIRE(SOMAGeometryDataFrame::exists(uri, ctx));
    REQUIRE(!SOMASparseNDArray::exists(uri, ctx));
    REQUIRE(!SOMADenseNDArray::exists(uri, ctx));
    REQUIRE(!SOMADataFrame::exists(uri, ctx));

    auto soma_geometry = SOMAGeometryDataFrame::open(
        uri, OpenMode::read, ctx, std::nullopt);
    REQUIRE(soma_geometry->uri() == uri);
    REQUIRE(soma_geometry->ctx() == ctx);
    REQUIRE(soma_geometry->type() == "SOMAGeometryDataFrame");
    std::vector<std::string> expected_index_column_names = {
        dim_infos[0].name, dim_infos[1].name};

    REQUIRE(soma_geometry->index_column_names() == expected_index_column_names);
    REQUIRE(soma_geometry->coordinate_space() == coord_space);
    REQUIRE(soma_geometry->nnz() == 0);
    soma_geometry->close();

    auto soma_object = SOMAObject::open(uri, OpenMode::read, ctx);
    REQUIRE(soma_object->uri() == uri);
    REQUIRE(soma_object->type() == "SOMAGeometryDataFrame");
    soma_object->close();
}

TEST_CASE("SOMAGeometryDataFrame: Roundtrip", "[SOMAGeometryDataFrame]") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri{"mem://unit-test-geometry-roundtrip"};
    PlatformConfig platform_config{};

    std::vector<helper::DimInfo> dim_infos(
        {helper::DimInfo(
             {.name = "soma_joinid",
              .tiledb_datatype = TILEDB_INT64,
              .dim_max = SOMA_JOINID_DIM_MAX,
              .string_lo = "N/A",
              .string_hi = "N/A"}),
         helper::DimInfo(
             {.name = "soma_geometry",
              .tiledb_datatype = TILEDB_GEOM_WKB,
              .dim_max = 100,
              .string_lo = "N/A",
              .string_hi = "N/A"})});

    std::vector<helper::AttrInfo> attr_infos({helper::AttrInfo(
        {.name = "quality", .tiledb_datatype = TILEDB_FLOAT64})});

    auto [schema, index_columns] =
        helper::create_arrow_schema_and_index_columns(
            dim_infos, attr_infos, coord_space);

    SOMAGeometryDataFrame::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        coord_space,
        ctx,
        platform_config,
        std::nullopt);

    // Create table of data for writing
    std::unique_ptr<ArrowSchema> data_schema = std::make_unique<ArrowSchema>(
        ArrowSchema{});
    std::unique_ptr<ArrowArray> data_array = std::make_unique<ArrowArray>(
        ArrowArray{});

    nanoarrow::UniqueBuffer metadata_buffer;
    ArrowMetadataBuilderInit(metadata_buffer.get(), nullptr);
    ArrowMetadataBuilderAppend(
        metadata_buffer.get(),
        ArrowCharView("geometry_type"),
        ArrowCharView("polygon_ring"));

    ArrowSchemaInitFromType(data_schema.get(), NANOARROW_TYPE_STRUCT);
    ArrowSchemaAllocateChildren(data_schema.get(), 3);
    ArrowSchemaInitFromType(data_schema->children[0], NANOARROW_TYPE_LIST);
    ArrowSchemaSetMetadata(
        data_schema->children[0],
        std::string((char*)metadata_buffer->data, metadata_buffer->size_bytes)
            .c_str());
    ArrowSchemaSetType(
        data_schema->children[0]->children[0], NANOARROW_TYPE_DOUBLE);
    ArrowSchemaSetName(data_schema->children[0], "soma_geometry");
    ArrowSchemaInitFromType(data_schema->children[1], NANOARROW_TYPE_INT64);
    ArrowSchemaSetName(data_schema->children[1], "soma_joinid");
    ArrowSchemaInitFromType(data_schema->children[2], NANOARROW_TYPE_DOUBLE);
    ArrowSchemaSetName(data_schema->children[2], "quality");

    ArrowArrayInitFromType(data_array.get(), NANOARROW_TYPE_STRUCT);
    ArrowArrayAllocateChildren(data_array.get(), 3);
    ArrowArrayInitFromType(data_array->children[0], NANOARROW_TYPE_LIST);
    ArrowArrayInitFromType(data_array->children[1], NANOARROW_TYPE_INT64);
    ArrowArrayInitFromType(data_array->children[2], NANOARROW_TYPE_DOUBLE);
    ArrowArrayAllocateChildren(data_array->children[0], 1);
    ArrowArrayInitFromType(
        data_array->children[0]->children[0], NANOARROW_TYPE_DOUBLE);
    ArrowArrayStartAppending(data_array->children[0]);
    ArrowArrayStartAppending(data_array->children[0]->children[0]);
    ArrowArrayStartAppending(data_array->children[1]);
    ArrowArrayStartAppending(data_array->children[2]);

    geometry::GenericGeometry polygon = geometry::Polygon(
        std::vector<geometry::BasePoint>(
            {geometry::BasePoint(0, 0),
             geometry::BasePoint(1, 0),
             geometry::BasePoint(0, 1)}));
    NANOARROW_THROW_NOT_OK(ArrowBufferAppendUInt32(
        ArrowArrayBuffer(data_array->children[0], 1), 0));
    data_array->children[0]->length = 1;
    NANOARROW_THROW_NOT_OK(
        ArrowArrayAppendDouble(data_array->children[0]->children[0], 0));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayAppendDouble(data_array->children[0]->children[0], 0));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayAppendDouble(data_array->children[0]->children[0], 1));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayAppendDouble(data_array->children[0]->children[0], 0));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayAppendDouble(data_array->children[0]->children[0], 0));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayAppendDouble(data_array->children[0]->children[0], 1));
    NANOARROW_THROW_NOT_OK(ArrowArrayAppendInt(data_array->children[1], 1));
    NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(data_array->children[2], 63));

    NANOARROW_THROW_NOT_OK(
        ArrowArrayFinishBuildingDefault(data_array->children[0], nullptr));
    NANOARROW_THROW_NOT_OK(ArrowArrayFinishBuildingDefault(
        data_array->children[0]->children[0], nullptr));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayFinishBuildingDefault(data_array->children[1], nullptr));
    NANOARROW_THROW_NOT_OK(
        ArrowArrayFinishBuildingDefault(data_array->children[2], nullptr));

    // Write to point cloud.
    auto soma_geometry = SOMAGeometryDataFrame::open(
        uri, OpenMode::write, ctx, std::nullopt);
    auto mq = ManagedQuery(*soma_geometry, ctx->tiledb_ctx());
    std::tie(data_array, data_schema) = TransformerPipeline(
                                            std::move(data_array),
                                            std::move(data_schema))
                                            .transform(
                                                OutlineTransformer(coord_space))
                                            .asTable();

    mq.set_array_data(std::move(data_schema), std::move(data_array));
    mq.submit_write();
    soma_geometry->close();

    // Read back the data.
    soma_geometry = SOMAGeometryDataFrame::open(
        uri, OpenMode::read, ctx, std::nullopt);
    mq = ManagedQuery(*soma_geometry, ctx->tiledb_ctx());
    while (auto batch = mq.read_next()) {
        auto arrbuf = batch.value();
        auto d0span = arrbuf->at(dim_infos[0].name)->data<int64_t>();
        auto d1span = arrbuf->at(SOMA_GEOMETRY_DIMENSION_PREFIX + "x__min")
                          ->data<double_t>();
        auto d2span = arrbuf->at(SOMA_GEOMETRY_DIMENSION_PREFIX + "x__max")
                          ->data<double_t>();
        auto d3span = arrbuf->at(SOMA_GEOMETRY_DIMENSION_PREFIX + "y__min")
                          ->data<double_t>();
        auto d4span = arrbuf->at(SOMA_GEOMETRY_DIMENSION_PREFIX + "y__max")
                          ->data<double_t>();
        auto wkbs = arrbuf->at(dim_infos[1].name)->binaries();
        auto a0span = arrbuf->at(attr_infos[0].name)->data<double>();
        CHECK(
            std::vector<int64_t>({1}) ==
            std::vector<int64_t>(d0span.begin(), d0span.end()));
        CHECK(
            std::vector<double_t>({0}) ==
            std::vector<double_t>(d1span.begin(), d1span.end()));
        CHECK(
            std::vector<double_t>({1}) ==
            std::vector<double_t>(d2span.begin(), d2span.end()));
        CHECK(
            std::vector<double_t>({0}) ==
            std::vector<double_t>(d3span.begin(), d3span.end()));
        CHECK(
            std::vector<double_t>({1}) ==
            std::vector<double_t>(d4span.begin(), d4span.end()));
        CHECK(geometry::to_wkb(polygon) == wkbs[0]);
        CHECK(
            std::vector<double_t>({63}) ==
            std::vector<double>(a0span.begin(), a0span.end()));
    }
    soma_geometry->close();

    auto soma_object = SOMAObject::open(uri, OpenMode::read, ctx);
    REQUIRE(soma_object->uri() == uri);
    REQUIRE(soma_object->type() == "SOMAGeometryDataFrame");
    soma_object->close();
}
