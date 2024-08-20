/**
 * @file   unit_soma_collection.cc
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
 * This file manages unit tests for the SOMACollection class
 */

#include "common.h"

static const int64_t DIM_MAX = 1000;

TEST_CASE("SOMACollection: basic") {
    TimestampRange ts(0, 2);
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-collection-basic";

    SOMACollection::create(uri, ctx, ts);
    auto soma_collection = SOMACollection::open(uri, OpenMode::read, ctx, ts);
    REQUIRE(soma_collection->uri() == uri);
    REQUIRE(soma_collection->ctx() == ctx);
    REQUIRE(soma_collection->type() == "SOMACollection");
    REQUIRE(soma_collection->timestamp() == ts);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMASparseNDArray") {
    bool use_current_domains[] = {false, true};
    for (bool use_current_domain : use_current_domains) {
        TimestampRange ts(0, 2);
        auto ctx = std::make_shared<SOMAContext>();
        std::string base_uri = "mem://unit-test-add-sparse-ndarray";
        std::string sub_uri = "mem://unit-test-add-sparse-ndarray/sub";

        SOMACollection::create(base_uri, ctx, ts);

        auto index_columns = helper::create_column_index_info(
            DIM_MAX, use_current_domain);

        std::map<std::string, SOMAGroupEntry> expected_map{
            {"sparse_ndarray", SOMAGroupEntry(sub_uri, "SOMAArray")}};

        auto soma_collection = SOMACollection::open(
            base_uri, OpenMode::write, ctx, ts);
        REQUIRE(soma_collection->timestamp() == ts);

        auto soma_sparse = soma_collection->add_new_sparse_ndarray(
            "sparse_ndarray",
            sub_uri,
            URIType::absolute,
            ctx,
            "l",
            ArrowTable(
                std::move(index_columns.first),
                std::move(index_columns.second)));
        REQUIRE(soma_collection->members_map() == expected_map);
        REQUIRE(soma_sparse->uri() == sub_uri);
        REQUIRE(soma_sparse->ctx() == ctx);
        REQUIRE(soma_sparse->type() == "SOMASparseNDArray");
        REQUIRE(soma_sparse->is_sparse() == true);
        REQUIRE(soma_sparse->ndim() == 1);
        REQUIRE(soma_sparse->nnz() == 0);
        REQUIRE(soma_sparse->timestamp() == ts);
        soma_sparse->close();
        soma_collection->close();

        soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
        REQUIRE(soma_collection->members_map() == expected_map);
        soma_collection->close();
    }
}

TEST_CASE("SOMACollection: add SOMADenseNDArray") {
    TimestampRange ts(0, 2);
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-add-dense-ndarray";
    std::string sub_uri = "mem://unit-test-add-dense-ndarray/sub";

    SOMACollection::create(base_uri, ctx, ts);
    auto index_columns = helper::create_column_index_info(DIM_MAX, false);

    std::map<std::string, SOMAGroupEntry> expected_map{
        {"dense_ndarray", SOMAGroupEntry(sub_uri, "SOMAArray")}};

    auto soma_collection = SOMACollection::open(
        base_uri, OpenMode::write, ctx, ts);
    REQUIRE(soma_collection->timestamp() == ts);

    auto soma_dense = soma_collection->add_new_dense_ndarray(
        "dense_ndarray",
        sub_uri,
        URIType::absolute,
        ctx,
        "l",
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)));
    REQUIRE(soma_collection->members_map() == expected_map);
    REQUIRE(soma_dense->uri() == sub_uri);
    REQUIRE(soma_dense->ctx() == ctx);
    REQUIRE(soma_dense->type() == "SOMADenseNDArray");
    REQUIRE(soma_dense->is_sparse() == false);
    REQUIRE(soma_dense->ndim() == 1);
    REQUIRE(soma_dense->shape() == std::vector<int64_t>{DIM_MAX + 1});
    REQUIRE(soma_dense->timestamp() == ts);
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->members_map() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMADataFrame") {
    TimestampRange ts(0, 2);
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-add-dataframe";
    std::string sub_uri = "mem://unit-test-add-dataframe/sub";

    SOMACollection::create(base_uri, ctx, ts);
    auto [schema, index_columns] = helper::create_arrow_schema(DIM_MAX, false);

    std::map<std::string, SOMAGroupEntry> expected_map{
        {"dataframe", SOMAGroupEntry(sub_uri, "SOMAArray")}};

    auto soma_collection = SOMACollection::open(
        base_uri, OpenMode::write, ctx, ts);
    REQUIRE(soma_collection->timestamp() == ts);

    auto soma_dataframe = soma_collection->add_new_dataframe(
        "dataframe",
        sub_uri,
        URIType::absolute,
        ctx,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)));
    REQUIRE(soma_collection->members_map() == expected_map);
    REQUIRE(soma_dataframe->uri() == sub_uri);
    REQUIRE(soma_dataframe->ctx() == ctx);
    REQUIRE(soma_dataframe->type() == "SOMADataFrame");
    std::vector<std::string> expected_index_column_names = {"d0"};
    REQUIRE(
        soma_dataframe->index_column_names() == expected_index_column_names);
    REQUIRE(soma_dataframe->timestamp() == ts);
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->members_map() == expected_map);
    REQUIRE(soma_dataframe->count() == 0);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMACollection") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-add-collection";
    std::string sub_uri = "mem://unit-test-add-collection/sub";

    SOMACollection::create(base_uri, ctx);

    std::map<std::string, SOMAGroupEntry> expected_map{
        {"subcollection", SOMAGroupEntry(sub_uri, "SOMAGroup")}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_subcollection = soma_collection->add_new_collection(
        "subcollection", sub_uri, URIType::absolute, ctx);
    REQUIRE(soma_collection->members_map() == expected_map);
    REQUIRE(soma_subcollection->uri() == sub_uri);
    REQUIRE(soma_subcollection->ctx() == ctx);
    REQUIRE(soma_subcollection->type() == "SOMACollection");
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->members_map() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMAExperiment") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-add-experiment";
    std::string sub_uri = "mem://unit-test-add-experiment/sub";

    SOMACollection::create(base_uri, ctx);
    auto [schema, index_columns] = helper::create_arrow_schema(DIM_MAX, false);

    std::map<std::string, SOMAGroupEntry> expected_map{
        {"experiment", SOMAGroupEntry(sub_uri, "SOMAGroup")}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_experiment = soma_collection->add_new_experiment(
        "experiment",
        sub_uri,
        URIType::absolute,
        ctx,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)));
    REQUIRE(soma_collection->members_map() == expected_map);
    REQUIRE(soma_experiment->uri() == sub_uri);
    REQUIRE(soma_experiment->ctx() == ctx);
    REQUIRE(soma_experiment->type() == "SOMAExperiment");
    soma_experiment->close();
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->members_map() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMAMeasurement") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-add-measurement";
    std::string sub_uri = "mem://unit-test-add-measurement/sub";

    SOMACollection::create(base_uri, ctx);
    auto [schema, index_columns] = helper::create_arrow_schema(DIM_MAX, false);

    std::map<std::string, SOMAGroupEntry> expected_map{
        {"measurement", SOMAGroupEntry(sub_uri, "SOMAGroup")}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_measurement = soma_collection->add_new_measurement(
        "measurement",
        sub_uri,
        URIType::absolute,
        ctx,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)));
    REQUIRE(soma_collection->members_map() == expected_map);
    REQUIRE(soma_measurement->uri() == sub_uri);
    REQUIRE(soma_measurement->ctx() == ctx);
    REQUIRE(soma_measurement->type() == "SOMAMeasurement");
    soma_measurement->close();
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->members_map() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: metadata") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-collection";
    SOMACollection::create(uri, ctx, TimestampRange(0, 2));
    auto soma_collection = SOMACollection::open(
        uri, OpenMode::write, ctx, std::pair<uint64_t, uint64_t>(1, 1));

    int32_t val = 100;
    soma_collection->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_collection->close();

    // Read metadata
    soma_collection->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(soma_collection->metadata_num() == 3);
    REQUIRE(soma_collection->has_metadata("soma_object_type"));
    REQUIRE(soma_collection->has_metadata("soma_encoding_version"));
    REQUIRE(soma_collection->has_metadata("md"));
    auto mdval = soma_collection->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_collection->close();

    // md should not be available at (2, 2)
    soma_collection->open(OpenMode::read, TimestampRange(2, 2));
    REQUIRE(soma_collection->metadata_num() == 2);
    REQUIRE(soma_collection->has_metadata("soma_object_type"));
    REQUIRE(soma_collection->has_metadata("soma_encoding_version"));
    REQUIRE(!soma_collection->has_metadata("md"));
    soma_collection->close();

    // Metadata should also be retrievable in write mode
    soma_collection->open(OpenMode::write, TimestampRange(0, 2));
    REQUIRE(soma_collection->metadata_num() == 3);
    REQUIRE(soma_collection->has_metadata("soma_object_type"));
    REQUIRE(soma_collection->has_metadata("soma_encoding_version"));
    REQUIRE(soma_collection->has_metadata("md"));
    mdval = soma_collection->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    soma_collection->delete_metadata("md");
    mdval = soma_collection->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_collection->close();

    // Confirm delete in read mode
    soma_collection->open(OpenMode::read, TimestampRange(0, 2));
    REQUIRE(!soma_collection->has_metadata("md"));
    REQUIRE(soma_collection->metadata_num() == 2);
}

TEST_CASE("SOMAExperiment: metadata") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-experiment";
    auto [schema, index_columns] = helper::create_arrow_schema(DIM_MAX, false);
    SOMAExperiment::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        PlatformConfig(),
        TimestampRange(0, 2));
    auto soma_experiment = SOMAExperiment::open(
        uri, OpenMode::write, ctx, std::pair<uint64_t, uint64_t>(1, 1));

    int32_t val = 100;
    soma_experiment->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_experiment->close();

    // Read metadata
    soma_experiment = SOMAExperiment::open(
        uri, OpenMode::read, ctx, TimestampRange(0, 2));
    REQUIRE(soma_experiment->metadata_num() == 4);
    REQUIRE(soma_experiment->has_metadata("dataset_type"));
    REQUIRE(soma_experiment->has_metadata("soma_object_type"));
    REQUIRE(soma_experiment->has_metadata("soma_encoding_version"));
    REQUIRE(soma_experiment->has_metadata("md"));
    auto mdval = soma_experiment->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_experiment->close();

    // md should not be available at (2, 2)
    soma_experiment = SOMAExperiment::open(
        uri, OpenMode::read, ctx, TimestampRange(2, 2));
    REQUIRE(soma_experiment->metadata_num() == 3);
    REQUIRE(soma_experiment->has_metadata("dataset_type"));
    REQUIRE(soma_experiment->has_metadata("soma_object_type"));
    REQUIRE(soma_experiment->has_metadata("soma_encoding_version"));
    REQUIRE(!soma_experiment->has_metadata("md"));
    soma_experiment->close();

    // Metadata should also be retrievable in write mode
    soma_experiment = SOMAExperiment::open(
        uri, OpenMode::write, ctx, TimestampRange(0, 2));
    REQUIRE(soma_experiment->metadata_num() == 4);
    REQUIRE(soma_experiment->has_metadata("dataset_type"));
    REQUIRE(soma_experiment->has_metadata("soma_object_type"));
    REQUIRE(soma_experiment->has_metadata("soma_encoding_version"));
    REQUIRE(soma_experiment->has_metadata("md"));
    mdval = soma_experiment->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    soma_experiment->delete_metadata("md");
    mdval = soma_experiment->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_experiment->close();

    // Confirm delete in read mode
    soma_experiment = SOMAExperiment::open(
        uri, OpenMode::read, ctx, TimestampRange(0, 2));
    REQUIRE(!soma_experiment->has_metadata("md"));
    REQUIRE(soma_experiment->metadata_num() == 3);
}

TEST_CASE("SOMAMeasurement: metadata") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string uri = "mem://unit-test-measurement";
    auto [schema, index_columns] = helper::create_arrow_schema(DIM_MAX, false);
    SOMAMeasurement::create(
        uri,
        std::move(schema),
        ArrowTable(
            std::move(index_columns.first), std::move(index_columns.second)),
        ctx,
        PlatformConfig(),
        TimestampRange(0, 2));

    auto soma_measurement = SOMAMeasurement::open(
        uri, OpenMode::write, ctx, std::pair<uint64_t, uint64_t>(1, 1));

    int32_t val = 100;
    soma_measurement->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_measurement->close();

    // Read metadata
    soma_measurement = SOMAMeasurement::open(
        uri, OpenMode::read, ctx, TimestampRange(0, 2));
    REQUIRE(soma_measurement->metadata_num() == 3);
    REQUIRE(soma_measurement->has_metadata("soma_object_type"));
    REQUIRE(soma_measurement->has_metadata("soma_encoding_version"));
    REQUIRE(soma_measurement->has_metadata("md"));
    auto mdval = soma_measurement->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_measurement->close();

    // md should not be available at (2, 2)
    soma_measurement = SOMAMeasurement::open(
        uri, OpenMode::read, ctx, TimestampRange(2, 2));
    REQUIRE(soma_measurement->metadata_num() == 2);
    REQUIRE(soma_measurement->has_metadata("soma_object_type"));
    REQUIRE(soma_measurement->has_metadata("soma_encoding_version"));
    REQUIRE(!soma_measurement->has_metadata("md"));
    soma_measurement->close();

    // Metadata should also be retrievable in write mode
    soma_measurement = SOMAMeasurement::open(
        uri, OpenMode::write, ctx, TimestampRange(0, 2));
    REQUIRE(soma_measurement->metadata_num() == 3);
    REQUIRE(soma_measurement->has_metadata("soma_object_type"));
    REQUIRE(soma_measurement->has_metadata("soma_encoding_version"));
    REQUIRE(soma_measurement->has_metadata("md"));
    mdval = soma_measurement->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);

    // Delete and have it reflected when reading metadata while in write mode
    soma_measurement->delete_metadata("md");
    mdval = soma_measurement->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_measurement->close();

    // Confirm delete in read mode
    soma_measurement = SOMAMeasurement::open(
        uri, OpenMode::read, ctx, TimestampRange(0, 2));
    REQUIRE(!soma_measurement->has_metadata("md"));
    REQUIRE(soma_measurement->metadata_num() == 2);
}
