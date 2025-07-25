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

#include "../utils/logger.h"
#include "soma_coordinates.h"
#include "tiledb_adapter/soma_query_condition.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMASparseNDArray::create(
    std::string_view uri,
    std::string_view format,
    const ArrowTable& index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    auto& index_column_schema = index_columns.second;
    uint64_t index_column_size = index_column_schema->n_children;

    auto schema = make_managed_unique<ArrowSchema>();
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

std::unique_ptr<SOMASparseNDArray> SOMASparseNDArray::open(
    std::string_view uri, OpenMode mode, std::shared_ptr<SOMAContext> ctx, std::optional<TimestampRange> timestamp) {
    auto array = std::make_unique<SOMASparseNDArray>(mode, uri, ctx, timestamp);

    if (!array->check_type("SOMASparseNDArray")) {
        throw TileDBSOMAError("[SOMASparseNDArray::open] Object is not a SOMASparseNDArray");
    }

    return array;
}

std::string_view SOMASparseNDArray::soma_data_type() {
    return ArrowAdapter::to_arrow_format(tiledb_schema()->attribute("soma_data").type());
}

//===================================================================
//= public non-static
//===================================================================

managed_unique_ptr<ArrowSchema> SOMASparseNDArray::schema() const {
    return this->arrow_schema();
}

void SOMASparseNDArray::delete_cells(const std::vector<SOMAColumnSelection<int64_t>>& coords) {
    if (coords.size() > ndim()) {
        throw std::invalid_argument(
            fmt::format(
                "Coordinates for {} columns were provided, but this array only has {} columns. The number of coords "
                "provided must be less than or equal to the number of columns.",
                coords.size(),
                ndim()));
    }
    const auto& array_shape = shape();

    SOMACoordQueryCondition qc{*ctx_, dimension_names()};
    for (size_t dim_index{0}; dim_index < coords.size(); ++dim_index) {
        qc.add_column_selection<int64_t>(
            dim_index, coords[dim_index], std::pair<int64_t, int64_t>(0, array_shape[dim_index] - 1));
    }

    auto soma_delete_cond = qc.get_soma_query_condition();
    if (!soma_delete_cond.is_initialized()) {
        throw std::invalid_argument("Cannot delete cells. At least one coordinate with values must be provided.");
    }
    delete_cells_impl(soma_delete_cond.query_condition());
}

}  // namespace tiledbsoma
