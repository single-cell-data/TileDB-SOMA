/**
 * @file   test_metadata.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the metadata layer
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

#include "common.h"

namespace {

std::tuple<std::string, uint64_t> create_array(
    const std::string& uri,
    std::shared_ptr<SOMAContext> ctx,
    int num_cells_per_fragment = 10,
    int num_fragments = 1,
    bool overlap = false,
    bool allow_duplicates = false) {
    auto vfs = VFS(*ctx->tiledb_ctx());
    if (vfs.is_dir(uri)) {
        vfs.remove_dir(uri);
    }

    const char* dim_name = "d0";
    const char* attr_name = "a0";

    // Create schema
    ArraySchema schema(*ctx->tiledb_ctx(), TILEDB_SPARSE);

    auto dim = Dimension::create<int64_t>(*ctx->tiledb_ctx(), dim_name, {0, std::numeric_limits<int64_t>::max() - 1});

    Domain domain(*ctx->tiledb_ctx());
    domain.add_dimension(dim);
    schema.set_domain(domain);

    auto attr = Attribute::create<int32_t>(*ctx->tiledb_ctx(), attr_name);
    schema.add_attribute(attr);
    schema.set_allows_dups(allow_duplicates);
    schema.check();

    // Create array
    SOMAArray::create(ctx, uri, std::move(schema), "NONE", std::nullopt);

    uint64_t nnz = num_fragments * num_cells_per_fragment;

    if (allow_duplicates) {
        return {uri, nnz};
    }

    // Adjust nnz when overlap is enabled
    if (overlap) {
        nnz = (num_fragments + 1) / 2 * num_cells_per_fragment;
    }

    return {uri, nnz};
}
}  // namespace

TEST_CASE("SOMAArray: metadata type check", "[Metadata]") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-array";
    const auto& [uri, expected_nnz] = create_array(base_uri, ctx);

    auto soma_array = SOMAArray::open(OpenMode::soma_write, uri, ctx);

    std::vector<int64_t> val_vec_int64({{-1, 100}});
    std::vector<uint64_t> val_vec_uint64({{1, 100}});
    std::vector<double> val_vec_double({{0.1, 100.01}});
    std::vector<int32_t> val_vec_int32({{-1, 100}});
    std::vector<uint32_t> val_vec_uint32({{1, 100}});
    std::vector<float> val_vec_float({{0.1, 100.01}});
    std::vector<int16_t> val_vec_int16({{-1, 100}});
    std::vector<uint16_t> val_vec_uint16({{1, 100}});
    std::vector<int8_t> val_vec_int8({{-1, 100}});
    std::vector<uint8_t> val_vec_uint8({{1, 100}});
    std::vector<uint8_t> val_vec_bool({{true, false}});
    std::vector<std::byte> val_vec_byte({{std::byte{1}, std::byte{100}}});
    std::string val_string = "test_string";
    int64_t val_int64 = -100;
    uint64_t val_uint64 = 100;
    double val_double = 100.01;
    int32_t val_int32 = -100;
    uint32_t val_uint32 = 100;
    float val_float = 100.01;
    int16_t val_int16 = -100;
    uint16_t val_uint16 = 100;
    int8_t val_int8 = -100;
    uint8_t val_uint8 = 100;
    uint8_t val_bool = true;

    soma_array->set_metadata("val_vec_int64", common::DataType::int64, val_vec_int64.size(), val_vec_int64.data());
    soma_array->set_metadata("val_vec_uint64", common::DataType::uint64, val_vec_uint64.size(), val_vec_uint64.data());
    soma_array->set_metadata("val_vec_double", common::DataType::float64, val_vec_double.size(), val_vec_double.data());
    soma_array->set_metadata("val_vec_int32", common::DataType::int32, val_vec_int32.size(), val_vec_int32.data());
    soma_array->set_metadata("val_vec_uint32", common::DataType::uint32, val_vec_uint32.size(), val_vec_uint32.data());
    soma_array->set_metadata("val_vec_float", common::DataType::float32, val_vec_float.size(), val_vec_float.data());
    soma_array->set_metadata("val_vec_int16", common::DataType::int16, val_vec_int16.size(), val_vec_int16.data());
    soma_array->set_metadata("val_vec_uint16", common::DataType::uint16, val_vec_uint16.size(), val_vec_uint16.data());
    soma_array->set_metadata("val_vec_int8", common::DataType::int8, val_vec_int8.size(), val_vec_int8.data());
    soma_array->set_metadata("val_vec_uint8", common::DataType::uint8, val_vec_uint8.size(), val_vec_uint8.data());
    // std::vector<bool> should not be used to write boolean data to TileDB due to implementation specific type size
    soma_array->set_metadata("val_vec_bool", common::DataType::boolean, val_vec_bool.size(), val_vec_bool.data());
    REQUIRE_THROWS(
        soma_array->set_metadata("val_vec_byte", common::DataType::blob, val_vec_byte.size(), val_vec_byte.data()));
    soma_array->set_metadata("val_string", common::DataType::string_utf8, val_string.size(), val_string.data());
    soma_array->set_metadata("val_int64", common::DataType::int64, 1, &val_int64);
    soma_array->set_metadata("val_uint64", common::DataType::uint64, 1, &val_uint64);
    soma_array->set_metadata("val_double", common::DataType::float64, 1, &val_double);
    soma_array->set_metadata("val_int32", common::DataType::int32, 1, &val_int32);
    soma_array->set_metadata("val_uint32", common::DataType::uint32, 1, &val_uint32);
    soma_array->set_metadata("val_float", common::DataType::float32, 1, &val_float);
    soma_array->set_metadata("val_int16", common::DataType::int16, 1, &val_int16);
    soma_array->set_metadata("val_uint16", common::DataType::uint16, 1, &val_uint16);
    soma_array->set_metadata("val_int8", common::DataType::int8, 1, &val_int8);
    soma_array->set_metadata("val_uint8", common::DataType::uint8, 1, &val_uint8);
    // bool should not be used to write boolean data to TileDB due to implementation specific type size
    soma_array->set_metadata("val_bool", common::DataType::boolean, 1, &val_bool);
    soma_array->close();

    // Read metadata
    soma_array->open(OpenMode::soma_read);
    REQUIRE(soma_array->metadata_num() == 25);
    REQUIRE(std::get<std::vector<int64_t>>(soma_array->get_metadata("val_vec_int64").value()) == val_vec_int64);
    REQUIRE(std::get<std::vector<uint64_t>>(soma_array->get_metadata("val_vec_uint64").value()) == val_vec_uint64);
    REQUIRE(std::get<std::vector<double>>(soma_array->get_metadata("val_vec_double").value()) == val_vec_double);
    REQUIRE(std::get<std::vector<int32_t>>(soma_array->get_metadata("val_vec_int32").value()) == val_vec_int32);
    REQUIRE(std::get<std::vector<uint32_t>>(soma_array->get_metadata("val_vec_uint32").value()) == val_vec_uint32);
    REQUIRE(std::get<std::vector<float>>(soma_array->get_metadata("val_vec_float").value()) == val_vec_float);
    REQUIRE(std::get<std::vector<int16_t>>(soma_array->get_metadata("val_vec_int16").value()) == val_vec_int16);
    REQUIRE(std::get<std::vector<uint16_t>>(soma_array->get_metadata("val_vec_uint16").value()) == val_vec_uint16);
    REQUIRE(std::get<std::vector<int8_t>>(soma_array->get_metadata("val_vec_int8").value()) == val_vec_int8);
    REQUIRE(std::get<std::vector<uint8_t>>(soma_array->get_metadata("val_vec_uint8").value()) == val_vec_uint8);
    // Reading boolean vector correctly transforms back to platform specific std::vector<bool>
    REQUIRE(
        std::get<std::vector<bool>>(soma_array->get_metadata("val_vec_bool").value()) ==
        std::vector<bool>(val_vec_bool.begin(), val_vec_bool.end()));
    REQUIRE(std::get<std::string>(soma_array->get_metadata("val_string").value()) == val_string);
    REQUIRE(std::get<int64_t>(soma_array->get_metadata("val_int64").value()) == val_int64);
    REQUIRE(std::get<uint64_t>(soma_array->get_metadata("val_uint64").value()) == val_uint64);
    REQUIRE(std::get<double>(soma_array->get_metadata("val_double").value()) == val_double);
    REQUIRE(std::get<int32_t>(soma_array->get_metadata("val_int32").value()) == val_int32);
    REQUIRE(std::get<uint32_t>(soma_array->get_metadata("val_uint32").value()) == val_uint32);
    REQUIRE(std::get<float>(soma_array->get_metadata("val_float").value()) == val_float);
    REQUIRE(std::get<int16_t>(soma_array->get_metadata("val_int16").value()) == val_int16);
    REQUIRE(std::get<uint16_t>(soma_array->get_metadata("val_uint16").value()) == val_uint16);
    REQUIRE(std::get<int8_t>(soma_array->get_metadata("val_int8").value()) == val_int8);
    REQUIRE(std::get<uint8_t>(soma_array->get_metadata("val_uint8").value()) == val_uint8);
    // Reading boolean vector correctly transforms back to platform specific std::vector<bool>
    REQUIRE(std::get<bool>(soma_array->get_metadata("val_bool").value()) == static_cast<bool>(val_bool));
    soma_array->close();
}

TEST_CASE("SOMAArray: metadata operations", "[Metadata]") {
    auto ctx = std::make_shared<SOMAContext>();
    std::string base_uri = "mem://unit-test-array";
    const auto& [uri, expected_nnz] = create_array(base_uri, ctx);

    auto soma_array = SOMAArray::open(OpenMode::soma_write, uri, ctx);

    int64_t val_wrong = 1000;
    int64_t val = -100;
    soma_array->set_metadata("val", common::DataType::int64, 1, &val_wrong);
    REQUIRE(soma_array->has_metadata("val"));
    soma_array->delete_metadata("val");
    REQUIRE(!soma_array->has_metadata("val"));
    soma_array->set_metadata("val", common::DataType::int64, 1, &val_wrong);
    soma_array->set_metadata("val", common::DataType::int64, 1, &val);
    REQUIRE(soma_array->has_metadata("val"));
    soma_array->close();

    // Read metadata
    soma_array->open(OpenMode::soma_read);
    REQUIRE(soma_array->metadata_num() == 3);
    REQUIRE(std::get<int64_t>(soma_array->get_metadata("val").value()) == val);
    REQUIRE_THROWS(soma_array->set_metadata("val", common::DataType::int64, 1, &val_wrong));
    REQUIRE_THROWS(soma_array->delete_metadata("val"));
    soma_array->close();
    // Metadata values are readable post close because the cache is not cleared
    REQUIRE(std::get<int64_t>(soma_array->get_metadata("val").value()) == val);
    REQUIRE_THROWS(soma_array->set_metadata("val", common::DataType::int64, 1, &val_wrong));
    REQUIRE_THROWS(soma_array->delete_metadata("val"));
}

TEST_CASE("SOMAGroup: metadata type check", "[Metadata]") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-group";
    SOMAGroup::create(ctx, uri, "NONE", {}, TimestampRange(0, 2));
    auto soma_group = SOMAGroup::open(OpenMode::soma_write, uri, ctx, "metadata", TimestampRange(1, 1));

    std::vector<int64_t> val_vec_int64({{-1, 100}});
    std::vector<uint64_t> val_vec_uint64({{1, 100}});
    std::vector<double> val_vec_double({{0.1, 100.01}});
    std::vector<int32_t> val_vec_int32({{-1, 100}});
    std::vector<uint32_t> val_vec_uint32({{1, 100}});
    std::vector<float> val_vec_float({{0.1, 100.01}});
    std::vector<int16_t> val_vec_int16({{-1, 100}});
    std::vector<uint16_t> val_vec_uint16({{1, 100}});
    std::vector<int8_t> val_vec_int8({{-1, 100}});
    std::vector<uint8_t> val_vec_uint8({{1, 100}});
    std::vector<uint8_t> val_vec_bool({{true, false}});
    std::vector<std::byte> val_vec_byte({{std::byte{1}, std::byte{100}}});
    std::string val_string = "test_string";
    int64_t val_int64 = -100;
    uint64_t val_uint64 = 100;
    double val_double = 100.01;
    int32_t val_int32 = -100;
    uint32_t val_uint32 = 100;
    float val_float = 100.01;
    int16_t val_int16 = -100;
    uint16_t val_uint16 = 100;
    int8_t val_int8 = -100;
    uint8_t val_uint8 = 100;
    uint8_t val_bool = true;

    soma_group->set_metadata("val_vec_int64", common::DataType::int64, val_vec_int64.size(), val_vec_int64.data());
    soma_group->set_metadata("val_vec_uint64", common::DataType::uint64, val_vec_uint64.size(), val_vec_uint64.data());
    soma_group->set_metadata("val_vec_double", common::DataType::float64, val_vec_double.size(), val_vec_double.data());
    soma_group->set_metadata("val_vec_int32", common::DataType::int32, val_vec_int32.size(), val_vec_int32.data());
    soma_group->set_metadata("val_vec_uint32", common::DataType::uint32, val_vec_uint32.size(), val_vec_uint32.data());
    soma_group->set_metadata("val_vec_float", common::DataType::float32, val_vec_float.size(), val_vec_float.data());
    soma_group->set_metadata("val_vec_int16", common::DataType::int16, val_vec_int16.size(), val_vec_int16.data());
    soma_group->set_metadata("val_vec_uint16", common::DataType::uint16, val_vec_uint16.size(), val_vec_uint16.data());
    soma_group->set_metadata("val_vec_int8", common::DataType::int8, val_vec_int8.size(), val_vec_int8.data());
    soma_group->set_metadata("val_vec_uint8", common::DataType::uint8, val_vec_uint8.size(), val_vec_uint8.data());
    // std::vector<bool> should not be used to write boolean data to TileDB due to implementation specific type size
    soma_group->set_metadata("val_vec_bool", common::DataType::boolean, val_vec_bool.size(), val_vec_bool.data());
    REQUIRE_THROWS(
        soma_group->set_metadata("val_vec_byte", common::DataType::blob, val_vec_byte.size(), val_vec_byte.data()));
    soma_group->set_metadata("val_string", common::DataType::string_utf8, val_string.size(), val_string.data());
    soma_group->set_metadata("val_int64", common::DataType::int64, 1, &val_int64);
    soma_group->set_metadata("val_uint64", common::DataType::uint64, 1, &val_uint64);
    soma_group->set_metadata("val_double", common::DataType::float64, 1, &val_double);
    soma_group->set_metadata("val_int32", common::DataType::int32, 1, &val_int32);
    soma_group->set_metadata("val_uint32", common::DataType::uint32, 1, &val_uint32);
    soma_group->set_metadata("val_float", common::DataType::float32, 1, &val_float);
    soma_group->set_metadata("val_int16", common::DataType::int16, 1, &val_int16);
    soma_group->set_metadata("val_uint16", common::DataType::uint16, 1, &val_uint16);
    soma_group->set_metadata("val_int8", common::DataType::int8, 1, &val_int8);
    soma_group->set_metadata("val_uint8", common::DataType::uint8, 1, &val_uint8);
    // bool should not be used to write boolean data to TileDB due to implementation specific type size
    soma_group->set_metadata("val_bool", common::DataType::boolean, 1, &val_bool);
    soma_group->close();

    // Read metadata
    soma_group->open(OpenMode::soma_read);
    REQUIRE(soma_group->metadata_num() == 23);
    REQUIRE(std::get<std::vector<int64_t>>(soma_group->get_metadata("val_vec_int64").value()) == val_vec_int64);
    REQUIRE(std::get<std::vector<uint64_t>>(soma_group->get_metadata("val_vec_uint64").value()) == val_vec_uint64);
    REQUIRE(std::get<std::vector<double>>(soma_group->get_metadata("val_vec_double").value()) == val_vec_double);
    REQUIRE(std::get<std::vector<int32_t>>(soma_group->get_metadata("val_vec_int32").value()) == val_vec_int32);
    REQUIRE(std::get<std::vector<uint32_t>>(soma_group->get_metadata("val_vec_uint32").value()) == val_vec_uint32);
    REQUIRE(std::get<std::vector<float>>(soma_group->get_metadata("val_vec_float").value()) == val_vec_float);
    REQUIRE(std::get<std::vector<int16_t>>(soma_group->get_metadata("val_vec_int16").value()) == val_vec_int16);
    REQUIRE(std::get<std::vector<uint16_t>>(soma_group->get_metadata("val_vec_uint16").value()) == val_vec_uint16);
    REQUIRE(std::get<std::vector<int8_t>>(soma_group->get_metadata("val_vec_int8").value()) == val_vec_int8);
    REQUIRE(std::get<std::vector<uint8_t>>(soma_group->get_metadata("val_vec_uint8").value()) == val_vec_uint8);
    // Reading boolean vector correctly transforms back to platform specific std::vector<bool>
    REQUIRE(
        std::get<std::vector<bool>>(soma_group->get_metadata("val_vec_bool").value()) ==
        std::vector<bool>(val_vec_bool.begin(), val_vec_bool.end()));
    REQUIRE(std::get<std::string>(soma_group->get_metadata("val_string").value()) == val_string);
    REQUIRE(std::get<int64_t>(soma_group->get_metadata("val_int64").value()) == val_int64);
    REQUIRE(std::get<uint64_t>(soma_group->get_metadata("val_uint64").value()) == val_uint64);
    REQUIRE(std::get<double>(soma_group->get_metadata("val_double").value()) == val_double);
    REQUIRE(std::get<int32_t>(soma_group->get_metadata("val_int32").value()) == val_int32);
    REQUIRE(std::get<uint32_t>(soma_group->get_metadata("val_uint32").value()) == val_uint32);
    REQUIRE(std::get<float>(soma_group->get_metadata("val_float").value()) == val_float);
    REQUIRE(std::get<int16_t>(soma_group->get_metadata("val_int16").value()) == val_int16);
    REQUIRE(std::get<uint16_t>(soma_group->get_metadata("val_uint16").value()) == val_uint16);
    REQUIRE(std::get<int8_t>(soma_group->get_metadata("val_int8").value()) == val_int8);
    REQUIRE(std::get<uint8_t>(soma_group->get_metadata("val_uint8").value()) == val_uint8);
    // Reading boolean vector correctly transforms back to platform specific std::vector<bool>
    REQUIRE(std::get<bool>(soma_group->get_metadata("val_bool").value()) == static_cast<bool>(val_bool));
    soma_group->close();
}

TEST_CASE("SOMAGroup: metadata operations", "[Metadata]") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-group";
    SOMAGroup::create(ctx, uri, "NONE", {}, TimestampRange(0, 2));
    auto soma_group = SOMAGroup::open(OpenMode::soma_write, uri, ctx, "metadata", TimestampRange(1, 1));

    int64_t val_wrong = 1000;
    int64_t val = -100;
    soma_group->set_metadata("val", common::DataType::int64, 1, &val_wrong);
    REQUIRE(soma_group->has_metadata("val"));
    soma_group->delete_metadata("val");
    REQUIRE(!soma_group->has_metadata("val"));
    soma_group->set_metadata("val", common::DataType::int64, 1, &val_wrong);
    soma_group->set_metadata("val", common::DataType::int64, 1, &val);
    REQUIRE(soma_group->has_metadata("val"));
    soma_group->close();

    // Read metadata
    soma_group->open(OpenMode::soma_read);
    REQUIRE(soma_group->metadata_num() == 1);
    REQUIRE(std::get<int64_t>(soma_group->get_metadata("val").value()) == val);
    REQUIRE_THROWS(soma_group->set_metadata("val", common::DataType::int64, 1, &val_wrong));
    REQUIRE_THROWS(soma_group->delete_metadata("val"));
    soma_group->close();
    // Metadata values are readable post close because the cache is not cleared
    REQUIRE(std::get<int64_t>(soma_group->get_metadata("val").value()) == val);
    REQUIRE_THROWS(soma_group->set_metadata("val", common::DataType::int64, 1, &val_wrong));
    REQUIRE_THROWS(soma_group->delete_metadata("val"));
}

TEST_CASE("Metadata edge cases", "[Metadata]") {
    auto ctx = std::make_shared<SOMAContext>();

    std::string uri = "mem://unit-test-group";
    SOMAGroup::create(ctx, uri, "NONE", {}, TimestampRange(0, 2));
    auto soma_group = SOMAGroup::open(OpenMode::soma_write, uri, ctx, "metadata", TimestampRange(1, 1));

    soma_group->set_metadata("empty string from python", common::DataType::string_utf8, 1, nullptr);
    soma_group->set_metadata("nullptr str", common::DataType::string_utf8, 0, nullptr);
    soma_group->set_metadata("nullptr int", common::DataType::int64, 0, nullptr);
    soma_group->set_metadata("nullptr bool", common::DataType::boolean, 0, nullptr);
    soma_group->close();

    // Read metadata
    soma_group->open(OpenMode::soma_read);
    REQUIRE(std::get<std::string>(soma_group->get_metadata("empty string from python").value()) == "");
    REQUIRE(std::get<std::string>(soma_group->get_metadata("nullptr str").value()) == "");
    REQUIRE(std::get<std::vector<int64_t>>(soma_group->get_metadata("nullptr int").value()) == std::vector<int64_t>());
    REQUIRE(std::get<std::vector<bool>>(soma_group->get_metadata("nullptr bool").value()) == std::vector<bool>());
    soma_group->close();
}