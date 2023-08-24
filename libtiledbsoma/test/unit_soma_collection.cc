/**
 * @file   unit_soma_collection.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_predicate.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include <numeric>
#include <random>

#include <tiledb/tiledb>
#include <tiledbsoma/tiledbsoma>
#include "utils/util.h"

using namespace tiledb;
using namespace tiledbsoma;
using namespace Catch::Matchers;

#ifndef TILEDBSOMA_SOURCE_ROOT
#define TILEDBSOMA_SOURCE_ROOT "not_defined"
#endif

const std::string src_path = TILEDBSOMA_SOURCE_ROOT;

namespace {
ArraySchema create_schema(
    Context& ctx, bool sparse = false, bool allow_duplicates = false) {
    // Create schema
    ArraySchema schema(ctx, sparse ? TILEDB_SPARSE : TILEDB_DENSE);

    auto dim = Dimension::create<int64_t>(ctx, "d0", {0, 1000});

    Domain domain(ctx);
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int>(ctx, "a0");
    schema.add_attribute(attr);
    schema.set_allows_dups(allow_duplicates);
    schema.check();

    return schema;
}
};  // namespace

TEST_CASE("SOMACollection: basic") {
    auto ctx = std::make_shared<Context>();
    std::string uri = "mem://unit-test-collection-basic";

    auto soma_collection = SOMACollection::create(uri, ctx);
    REQUIRE(soma_collection->uri() == uri);
    REQUIRE(soma_collection->ctx() == ctx);
    REQUIRE(soma_collection->type() == "SOMACollection");
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMASparseNDArray") {
    auto ctx = std::make_shared<Context>();
    std::string base_uri = "mem://unit-test-add-sparse-ndarray";
    std::string sub_uri = "mem://unit-test-add-sparse-ndarray/sub";

    SOMACollection::create(base_uri, ctx);
    auto schema = create_schema(*ctx, true);

    std::map<std::string, std::string> expected_map{
        {"sparse_ndarray", sub_uri}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_sparse = soma_collection->add_new_sparse_ndarray(
        "sparse_ndarray", sub_uri, URIType::absolute, ctx, schema);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    REQUIRE(soma_sparse->uri() == sub_uri);
    REQUIRE(soma_sparse->ctx() == ctx);
    REQUIRE(soma_sparse->type() == "SOMASparseNDArray");
    REQUIRE(soma_sparse->is_sparse() == true);
    REQUIRE(soma_sparse->schema()->has_attribute("a0"));
    REQUIRE(soma_sparse->schema()->domain().has_dimension("d0"));
    REQUIRE(soma_sparse->ndim() == 1);
    REQUIRE(soma_sparse->nnz() == 0);
    soma_sparse->close();
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMADenseNDArray") {
    auto ctx = std::make_shared<Context>();
    std::string base_uri = "mem://unit-test-add-dense-ndarray";
    std::string sub_uri = "mem://unit-test-add-dense-ndarray/sub";

    SOMACollection::create(base_uri, ctx);
    auto schema = create_schema(*ctx, false);

    std::map<std::string, std::string> expected_map{{"dense_ndarray", sub_uri}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_dense = soma_collection->add_new_dense_ndarray(
        "dense_ndarray", sub_uri, URIType::absolute, ctx, schema);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    REQUIRE(soma_dense->uri() == sub_uri);
    REQUIRE(soma_dense->ctx() == ctx);
    REQUIRE(soma_dense->type() == "SOMADenseNDArray");
    REQUIRE(soma_dense->is_sparse() == false);
    REQUIRE(soma_dense->schema()->has_attribute("a0"));
    REQUIRE(soma_dense->schema()->domain().has_dimension("d0"));
    REQUIRE(soma_dense->ndim() == 1);
    REQUIRE(soma_dense->shape() == std::vector<int64_t>{1001});
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMADataFrame") {
    auto ctx = std::make_shared<Context>();
    std::string base_uri = "mem://unit-test-add-dataframe";
    std::string sub_uri = "mem://unit-test-add-dataframe/sub";

    SOMACollection::create(base_uri, ctx);
    auto schema = create_schema(*ctx, false);

    std::map<std::string, std::string> expected_map{{"dataframe", sub_uri}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_dataframe = soma_collection->add_new_dataframe(
        "dataframe", sub_uri, URIType::absolute, ctx, schema);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    REQUIRE(soma_dataframe->uri() == sub_uri);
    REQUIRE(soma_dataframe->ctx() == ctx);
    REQUIRE(soma_dataframe->type() == "SOMADataFrame");
    REQUIRE(soma_dataframe->schema()->has_attribute("a0"));
    REQUIRE(soma_dataframe->schema()->domain().has_dimension("d0"));
    std::vector<std::string> expected_index_column_names = {"d0"};
    REQUIRE(
        soma_dataframe->index_column_names() == expected_index_column_names);
    REQUIRE(soma_dataframe->count() == 1);
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMACollection") {
    auto ctx = std::make_shared<Context>();
    std::string base_uri = "mem://unit-test-add-collection";
    std::string sub_uri = "mem://unit-test-add-collection/sub";

    SOMACollection::create(base_uri, ctx);
    auto schema = create_schema(*ctx, false);

    std::map<std::string, std::string> expected_map{{"subcollection", sub_uri}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_subcollection = soma_collection->add_new_collection(
        "subcollection", sub_uri, URIType::absolute, ctx);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    REQUIRE(soma_subcollection->uri() == sub_uri);
    REQUIRE(soma_subcollection->ctx() == ctx);
    REQUIRE(soma_subcollection->type() == "SOMACollection");
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMAExperiment") {
    auto ctx = std::make_shared<Context>();
    std::string base_uri = "mem://unit-test-add-experiment";
    std::string sub_uri = "mem://unit-test-add-experiment/sub";

    SOMACollection::create(base_uri, ctx);
    auto schema = create_schema(*ctx, false);

    std::map<std::string, std::string> expected_map{{"experiment", sub_uri}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_experiment = soma_collection->add_new_experiment(
        "experiment", sub_uri, URIType::absolute, ctx, schema);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    REQUIRE(soma_experiment->uri() == sub_uri);
    REQUIRE(soma_experiment->ctx() == ctx);
    REQUIRE(soma_experiment->type() == "SOMAExperiment");
    soma_experiment->close();
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: add SOMAMeasurement") {
    auto ctx = std::make_shared<Context>();
    std::string base_uri = "mem://unit-test-add-measurement";
    std::string sub_uri = "mem://unit-test-add-measurement/sub";

    SOMACollection::create(base_uri, ctx);
    auto schema = create_schema(*ctx, false);

    std::map<std::string, std::string> expected_map{{"measurement", sub_uri}};

    auto soma_collection = SOMACollection::open(base_uri, OpenMode::write, ctx);
    auto soma_measurement = soma_collection->add_new_measurement(
        "measurement", sub_uri, URIType::absolute, ctx, schema);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    REQUIRE(soma_measurement->uri() == sub_uri);
    REQUIRE(soma_measurement->ctx() == ctx);
    REQUIRE(soma_measurement->type() == "SOMAMeasurement");
    soma_measurement->close();
    soma_collection->close();

    soma_collection = SOMACollection::open(base_uri, OpenMode::read, ctx);
    REQUIRE(soma_collection->member_to_uri_mapping() == expected_map);
    soma_collection->close();
}

TEST_CASE("SOMACollection: metadata") {
    auto ctx = std::make_shared<Context>();

    std::string uri = "mem://unit-test-collection";
    SOMACollection::create(uri, ctx);
    auto soma_collection = SOMACollection::open(
        uri, OpenMode::write, ctx, std::pair<uint64_t, uint64_t>(1, 1));
    int32_t val = 100;
    soma_collection->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_collection->close();

    soma_collection->open(OpenMode::read, std::pair<uint64_t, uint64_t>(1, 1));
    REQUIRE(soma_collection->metadata_num() == 2);
    REQUIRE(soma_collection->has_metadata("soma_object_type") == true);
    REQUIRE(soma_collection->has_metadata("md") == true);

    auto mdval = soma_collection->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_collection->close();

    soma_collection->open(OpenMode::write, std::pair<uint64_t, uint64_t>(2, 2));
    // Metadata should also be retrievable in write mode
    mdval = soma_collection->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_collection->delete_metadata("md");
    mdval = soma_collection->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_collection->close();

    soma_collection->open(OpenMode::read, std::pair<uint64_t, uint64_t>(3, 3));
    REQUIRE(soma_collection->has_metadata("md") == false);
    REQUIRE(soma_collection->metadata_num() == 1);
    soma_collection->close();
}

TEST_CASE("SOMAExperiment: metadata") {
    auto ctx = std::make_shared<Context>();

    std::string uri = "mem://unit-test-experiment";
    SOMAExperiment::create(uri, create_schema(*ctx), ctx);
    auto soma_experiment = SOMAExperiment::open(
        uri, OpenMode::write, ctx, std::pair<uint64_t, uint64_t>(1, 1));
    int32_t val = 100;
    soma_experiment->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_experiment->close();

    soma_experiment->open(OpenMode::read, std::pair<uint64_t, uint64_t>(1, 1));
    REQUIRE(soma_experiment->metadata_num() == 2);
    REQUIRE(soma_experiment->has_metadata("soma_object_type") == true);
    REQUIRE(soma_experiment->has_metadata("md") == true);

    auto mdval = soma_experiment->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_experiment->close();

    soma_experiment->open(OpenMode::write, std::pair<uint64_t, uint64_t>(2, 2));
    // Metadata should also be retrievable in write mode
    mdval = soma_experiment->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_experiment->delete_metadata("md");
    mdval = soma_experiment->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_experiment->close();

    soma_experiment->open(OpenMode::read, std::pair<uint64_t, uint64_t>(3, 3));
    REQUIRE(soma_experiment->has_metadata("md") == false);
    REQUIRE(soma_experiment->metadata_num() == 1);
    soma_experiment->close();
}

TEST_CASE("SOMAMeasurement: metadata") {
    auto ctx = std::make_shared<Context>();

    std::string uri = "mem://unit-test-measurement";
    SOMAMeasurement::create(uri, create_schema(*ctx), ctx);
    auto soma_measurement = SOMAMeasurement::open(
        uri, OpenMode::write, ctx, std::pair<uint64_t, uint64_t>(1, 1));
    int32_t val = 100;
    soma_measurement->set_metadata("md", TILEDB_INT32, 1, &val);
    soma_measurement->close();

    soma_measurement->open(OpenMode::read, std::pair<uint64_t, uint64_t>(1, 1));
    REQUIRE(soma_measurement->metadata_num() == 2);
    REQUIRE(soma_measurement->has_metadata("soma_object_type") == true);
    REQUIRE(soma_measurement->has_metadata("md") == true);

    auto mdval = soma_measurement->get_metadata("md");
    REQUIRE(std::get<MetadataInfo::dtype>(*mdval) == TILEDB_INT32);
    REQUIRE(std::get<MetadataInfo::num>(*mdval) == 1);
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_measurement->close();

    soma_measurement->open(
        OpenMode::write, std::pair<uint64_t, uint64_t>(2, 2));
    // Metadata should also be retrievable in write mode
    mdval = soma_measurement->get_metadata("md");
    REQUIRE(*((const int32_t*)std::get<MetadataInfo::value>(*mdval)) == 100);
    soma_measurement->delete_metadata("md");
    mdval = soma_measurement->get_metadata("md");
    REQUIRE(!mdval.has_value());
    soma_measurement->close();

    soma_measurement->open(OpenMode::read, std::pair<uint64_t, uint64_t>(3, 3));
    REQUIRE(soma_measurement->has_metadata("md") == false);
    REQUIRE(soma_measurement->metadata_num() == 1);
    soma_measurement->close();
}