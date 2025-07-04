/**
 * @file   unit_arrow_adapter.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages unit tests for the SOMADataFrame class
 */

#include "common.h"

TEST_CASE("name", "[pattern]") {
    std::vector<std::string> names({"int64", "uint8", "float64", "float32", "bool", "string"});
    std::vector<tiledb_datatype_t> tiledb_datatypes(
        {TILEDB_INT64, TILEDB_UINT8, TILEDB_FLOAT64, TILEDB_FLOAT32, TILEDB_BOOL, TILEDB_STRING_ASCII});

    REQUIRE(names.size() == tiledb_datatypes.size());
    auto n_columns = names.size();

    managed_unique_ptr<ArrowSchema> arrow_schema = ArrowAdapter::make_arrow_schema(names, tiledb_datatypes);
    managed_unique_ptr<ArrowArray> arrow_array = ArrowAdapter::make_arrow_array_parent(n_columns);

    std::vector<int64_t> inputs_int64({1000, 2000, 3000});
    std::vector<uint8_t> inputs_uint8({33, 44, 55});
    std::vector<double> inputs_float64({1.0, 2.5, 4.0});
    std::vector<float> inputs_float32({10.0, 12.5, 14.0});
    std::vector<bool> inputs_bool({true, false, true});
    std::vector<std::string> inputs_string({"apple", "", "cat"});

    arrow_array->children[0] = ArrowAdapter::make_arrow_array_child(inputs_int64);
    arrow_array->children[1] = ArrowAdapter::make_arrow_array_child(inputs_uint8);
    arrow_array->children[2] = ArrowAdapter::make_arrow_array_child(inputs_float64);
    arrow_array->children[3] = ArrowAdapter::make_arrow_array_child(inputs_float32);
    arrow_array->children[4] = ArrowAdapter::make_arrow_array_child(inputs_bool);
    arrow_array->children[5] = ArrowAdapter::make_arrow_array_child_string(inputs_string);

    ArrowTable arrow_table(std::move(arrow_array), std::move(arrow_schema));

    std::vector<int64_t> outputs_int64 = ArrowAdapter::get_table_non_string_column_by_name<int64_t>(
        arrow_table, "int64");
    std::vector<uint8_t> outputs_uint8 = ArrowAdapter::get_table_non_string_column_by_name<uint8_t>(
        arrow_table, "uint8");
    std::vector<double> outputs_float64 = ArrowAdapter::get_table_non_string_column_by_name<double>(
        arrow_table, "float64");
    std::vector<float> outputs_float32 = ArrowAdapter::get_table_non_string_column_by_name<float>(
        arrow_table, "float32");
    std::vector<bool> outputs_bool = ArrowAdapter::get_table_non_string_column_by_name<bool>(arrow_table, "bool");
    std::vector<std::string> outputs_string = ArrowAdapter::get_table_string_column_by_name(arrow_table, "string");

    REQUIRE(outputs_int64 == inputs_int64);
    REQUIRE(outputs_uint8 == inputs_uint8);
    REQUIRE(outputs_float64 == inputs_float64);
    REQUIRE(outputs_float32 == inputs_float32);
    REQUIRE(outputs_bool == inputs_bool);
    REQUIRE(outputs_string == inputs_string);
}
