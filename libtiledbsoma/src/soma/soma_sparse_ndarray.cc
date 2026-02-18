/**
 * @file   soma_sparse_ndarray.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMASparseNDArray class.
 */

#include "soma_sparse_ndarray.h"

#include <stdexcept>
#include <type_traits>

#include "../utils/arrow_adapter.h"
#include "common/logging/impl/logger.h"
#include "soma_coordinates.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMASparseNDArray::create(
    std::string_view uri,
    std::string_view format,
    const common::arrow::ArrowTable& index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    auto& index_column_schema = index_columns.second;
    uint64_t index_column_size = index_column_schema->n_children;

    auto schema = common::arrow::make_managed_unique<ArrowSchema>();
    schema->name = nullptr;
    schema->format = strdup("+s");
    schema->n_children = index_column_size + 1;
    schema->dictionary = nullptr;
    schema->metadata = nullptr;
    schema->flags = 0;
    schema->release = &ArrowAdapter::release_schema;
    schema->children = (ArrowSchema**)malloc(schema->n_children * sizeof(ArrowSchema*));

    std::vector<std::string> index_column_names;
    for (uint64_t dim_idx = 0; dim_idx < index_column_size; ++dim_idx) {
        ArrowSchema* dim = schema->children[dim_idx] = (ArrowSchema*)malloc(sizeof(ArrowSchema));
        dim->format = strdup("l");
        dim->name = strdup(std::string("soma_dim_" + std::to_string(dim_idx)).c_str());
        dim->n_children = 0;
        dim->children = nullptr;
        dim->dictionary = nullptr;
        dim->metadata = nullptr;
        dim->release = &ArrowAdapter::release_schema;
        index_column_names.push_back(dim->name);
    }

    ArrowSchema* attr = schema->children[index_column_size] = (ArrowSchema*)malloc(sizeof(ArrowSchema));
    attr->format = strdup(std::string(format).c_str());
    attr->name = strdup("soma_data");
    attr->n_children = 0;
    attr->flags = 0;
    attr->children = nullptr;
    attr->dictionary = nullptr;
    attr->metadata = nullptr;
    attr->release = &ArrowAdapter::release_schema;

    auto [tiledb_schema, soma_schema_extension] = ArrowAdapter::tiledb_schema_from_arrow_schema(
        ctx->tiledb_ctx(), schema, index_columns, std::nullopt, "SOMASparseNDArray", true, platform_config, timestamp);

    SOMAArray::create(ctx, uri, tiledb_schema, "SOMASparseNDArray", std::nullopt, timestamp);
}

void SOMASparseNDArray::create(
    std::string_view uri,
    std::string_view format,
    std::span<const std::optional<int64_t>> shape,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    std::vector<int64_t> sanitized_shape;
    std::transform(shape.begin(), shape.end(), std::back_inserter(sanitized_shape), [](const auto& dim_shape) {
        return dim_shape.value_or(1);
    });

    tiledb::ArraySchema schema = utils::create_nd_array_schema(
        "SOMASparseNDArray", true, format, sanitized_shape, ctx->tiledb_ctx(), platform_config, timestamp);

    SOMAArray::create(ctx, uri, schema, "SOMASparseNDArray", std::nullopt, timestamp);
}

std::unique_ptr<SOMASparseNDArray> SOMASparseNDArray::open(
    std::string_view uri, OpenMode mode, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    return std::make_unique<SOMASparseNDArray>(mode, uri, ctx, timestamp);
}

std::string_view SOMASparseNDArray::soma_data_type() {
    return ArrowAdapter::to_arrow_format(tiledb_schema()->attribute("soma_data").type());
}

//===================================================================
//= public non-static
//===================================================================

common::arrow::managed_unique_ptr<ArrowSchema> SOMASparseNDArray::schema() const {
    return this->arrow_schema(true);
}

}  // namespace tiledbsoma
