/**
 * @file   soma_dense_ndarray.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMADenseNDArray class.
 */
#include "soma_dense_ndarray.h"
#include "soma_coordinates.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

void SOMADenseNDArray::create(
    std::string_view uri,
    std::string_view format,
    const ArrowTable& index_columns,
    std::shared_ptr<SOMAContext> ctx,
    PlatformConfig platform_config,
    std::optional<TimestampRange> timestamp) {
    auto& index_column_schema = index_columns.second;
    uint64_t index_column_size = index_column_schema->n_children;

    auto schema = std::make_unique<ArrowSchema>();
    schema->format = strdup("+s");
    schema->n_children = index_column_size + 1;
    schema->dictionary = nullptr;
    schema->metadata = nullptr;
    schema->flags = 0;
    schema->release = &ArrowAdapter::release_schema;
    schema->children = new ArrowSchema*[schema->n_children];

    std::vector<std::string> index_column_names;
    for (uint64_t dim_idx = 0; dim_idx < index_column_size; ++dim_idx) {
        ArrowSchema* dim = schema->children[dim_idx] = new ArrowSchema;
        dim->format = strdup("l");
        dim->name = strdup(
            std::string("soma_dim_" + std::to_string(dim_idx)).c_str());
        dim->n_children = 0;
        dim->metadata = nullptr;
        dim->dictionary = nullptr;
        dim->release = &ArrowAdapter::release_schema;
        index_column_names.push_back(dim->name);
    }

    ArrowSchema* attr = schema->children[index_column_size] = new ArrowSchema;
    attr->format = strdup(std::string(format).c_str());
    attr->name = strdup("soma_data");
    attr->n_children = 0;
    attr->flags = 0;  // or ARROW_FLAG_NULLABLE;
    attr->dictionary = nullptr;
    attr->metadata = nullptr;
    attr->release = &ArrowAdapter::release_schema;

    auto [tiledb_schema, soma_schema_extension] =
        ArrowAdapter::tiledb_schema_from_arrow_schema(
            ctx->tiledb_ctx(),
            schema,
            index_columns,
            std::nullopt,
            "SOMADenseNDArray",
            false,
            platform_config);

    SOMAArray::create(
        ctx, uri, tiledb_schema, "SOMADenseNDArray", std::nullopt, timestamp);
}

std::unique_ptr<SOMADenseNDArray> SOMADenseNDArray::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::vector<std::string> column_names,
    ResultOrder result_order,
    std::optional<TimestampRange> timestamp) {
    auto array = std::make_unique<SOMADenseNDArray>(
        mode, uri, ctx, column_names, result_order, timestamp);

    if (!array->check_type("SOMADenseNDArray")) {
        throw TileDBSOMAError(
            "[SOMADenseNDArray::open] Object is not a SOMADenseNDArray");
    }

    return array;
}

bool SOMADenseNDArray::exists(
    std::string_view uri, std::shared_ptr<SOMAContext> ctx) {
    try {
        auto obj = SOMAObject::open(uri, OpenMode::read, ctx);
        return "SOMADenseNDArray" == obj->type();
    } catch (TileDBSOMAError& e) {
        return false;
    }
}

std::string_view SOMADenseNDArray::soma_data_type() {
    return ArrowAdapter::to_arrow_format(
        tiledb_schema()->attribute("soma_data").type());
}

//===================================================================
//= public non-static
//===================================================================

std::unique_ptr<ArrowSchema> SOMADenseNDArray::schema() const {
    return this->arrow_schema();
}

}  // namespace tiledbsoma
