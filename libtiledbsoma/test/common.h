/**
 * @file   common.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages common headers and helper classes for the unit test files.
 */

#ifndef UNIT_TEST_COMMON_H
#define UNIT_TEST_COMMON_H

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

static const std::string src_path = TILEDBSOMA_SOURCE_ROOT;

namespace helper {

// This non-obvious number is:
// * Something that fits into signed 32-bit integer for R-friendliness;
// * Is a comfortable tile-extent distance away from 2^31-1 for default
//   core tile extent. (Using 2^31-1 exactly would result in a core
//   array-creation error.)
const int CORE_DOMAIN_MAX = 2147483646;

// E.g. "d0" is of type TILEDB_INT64 with dim_max 1000 and current-domain
// feature enabled
struct DimInfo {
    std::string name;
    tiledb_datatype_t tiledb_datatype;
    int64_t dim_max;
    std::string string_lo;  // For custom/restricted DataFrame domains
    std::string string_hi;  // For custom/restricted DataFrame domains
};

// E.g. "a0" is of type TILEDB_FLOAT64
struct AttrInfo {
    std::string name;
    tiledb_datatype_t tiledb_datatype;
};

std::pair<std::unique_ptr<ArrowSchema>, ArrowTable>
create_arrow_schema_and_index_columns(
    const std::vector<DimInfo>& dim_infos,
    const std::vector<AttrInfo>& attr_infos,
    std::optional<SOMACoordinateSpace> coordinate_space = std::nullopt);

ArrowTable create_column_index_info(const std::vector<DimInfo>& dim_infos);

std::string to_arrow_format(tiledb_datatype_t tiledb_datatype);

// Core PR: https://github.com/TileDB-Inc/TileDB/pull/5303
bool have_dense_current_domain_support();

std::unique_ptr<ArrowSchema> create_index_cols_info_schema(
    const std::vector<DimInfo>& dim_infos,
    std::optional<SOMACoordinateSpace> coordinate_space = std::nullopt);

}  // namespace helper
#endif
