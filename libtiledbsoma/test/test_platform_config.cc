/**
 * @file  test_soma_column_selection.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file manages tests for the SOMA column selection classes.
 */

#include "common.h"

TEMPLATE_TEST_CASE_SIG(
    "DimensionConfigAdapter: basic signed type with default tiledb type",
    "[DimensionConfigAdpater][PlatformConfig]",
    ((typename T, tiledb_datatype_t D), T, D),
    (int8_t, TILEDB_INT8),
    (int16_t, TILEDB_INT16),
    (int32_t, TILEDB_INT32),
    (int64_t, TILEDB_INT64),
    (float_t, TILEDB_FLOAT32),
    (double_t, TILEDB_FLOAT64)) {
    std::array<T, 2> max_domain{10, 100};
    std::array<T, 2> current_domain{10, 50};
    T tile_extent{20};
    DimensionConfigAdapter<T> dim_config(current_domain, max_domain, tile_extent);

    Context ctx{};
    std::string name{"index_column_0"};
    FilterList filter_list{ctx};
    tiledb::Dimension dim = dim_config.create_dimension(ctx, name, filter_list);
    CHECK(dim.name() == name);
    REQUIRE(dim.type() == D);
    auto actual_max_domain = dim.domain<T>();
    auto actual_tile_extent = dim.tile_extent<T>();
    CHECK(actual_max_domain.first == max_domain[0]);
    CHECK(actual_max_domain.second == max_domain[1]);
    CHECK(actual_tile_extent == tile_extent);
}

TEMPLATE_TEST_CASE_SIG(
    "DimensionConfigAdapter: basic signed types with specific tiledb types",
    "[DimensionConfigAdpater][PlatformConfig]",
    ((typename T, tiledb_datatype_t D), T, D),
    (int8_t, TILEDB_INT8),
    (int16_t, TILEDB_INT16),
    (int32_t, TILEDB_INT32),
    (int64_t, TILEDB_INT64),
    (int64_t, TILEDB_DATETIME_SEC),
    (int64_t, TILEDB_DATETIME_MS),
    (int64_t, TILEDB_DATETIME_US),
    (int64_t, TILEDB_DATETIME_NS),
    (float_t, TILEDB_FLOAT32),
    (double_t, TILEDB_FLOAT64)) {
    std::array<T, 2> max_domain{10, 100};
    std::array<T, 2> current_domain{10, 50};
    T tile_extent{20};
    DimensionConfigAdapter<T> dim_config(current_domain, max_domain, tile_extent);

    Context ctx{};
    std::string name{"index_column_0"};
    FilterList filter_list{ctx};
    tiledb::Dimension dim = dim_config.create_dimension(ctx, name, D, filter_list);
    CHECK(dim.name() == name);
    auto actual_max_domain = dim.domain<T>();
    auto actual_tile_extent = dim.tile_extent<T>();
    CHECK(actual_max_domain.first == max_domain[0]);
    CHECK(actual_max_domain.second == max_domain[1]);
    CHECK(actual_tile_extent == tile_extent);
}
