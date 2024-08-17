/**
 * @file   common.cc
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
 * This file manages common headers and helper classes for the unit test files.
 */

#include "common.h"

namespace helper {
ArraySchema create_schema(Context& ctx, bool allow_duplicates) {
    // Create schema
    ArraySchema schema(ctx, TILEDB_SPARSE);

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

std::pair<std::unique_ptr<ArrowSchema>, ArrowTable> create_arrow_schema() {
    // Create ArrowSchema for SOMAArray
    auto arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = "+s";
    arrow_schema->n_children = 2;
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->children = new ArrowSchema*[arrow_schema->n_children];
    ArrowSchema* dim = arrow_schema->children[0] = new ArrowSchema;
    dim->format = "l";
    dim->name = "d0";
    dim->n_children = 0;
    dim->dictionary = nullptr;
    dim->release = &ArrowAdapter::release_schema;
    ArrowSchema* attr = arrow_schema->children[1] = new ArrowSchema;
    attr->format = "l";
    attr->name = "a0";
    attr->n_children = 0;
    attr->flags = 0;
    attr->dictionary = nullptr;
    attr->release = &ArrowAdapter::release_schema;

    // Create ArrowSchema for IndexColumnInfo
    auto col_info_schema = std::make_unique<ArrowSchema>();
    col_info_schema->format = "+s";
    col_info_schema->n_children = 1;
    col_info_schema->dictionary = nullptr;
    col_info_schema->release = &ArrowAdapter::release_schema;
    col_info_schema->children = new ArrowSchema*[col_info_schema->n_children];
    dim = col_info_schema->children[0] = new ArrowSchema;
    dim->format = "l";
    dim->name = "d0";
    dim->n_children = 0;
    dim->dictionary = nullptr;
    dim->release = &ArrowAdapter::release_schema;

    // Create ArrowArray for IndexColumnInfo
    auto col_info_array = std::make_unique<ArrowArray>();
    col_info_array->length = 0;
    col_info_array->null_count = 0;
    col_info_array->offset = 0;
    col_info_array->n_buffers = 0;
    col_info_array->buffers = nullptr;
    col_info_array->n_children = 2;
    col_info_array->release = &ArrowAdapter::release_array;
    col_info_array->children = new ArrowArray*[1];
    int n = 3;
    auto d0_info = col_info_array->children[0] = new ArrowArray;
    d0_info->length = n;
    d0_info->null_count = 0;
    d0_info->offset = 0;
    d0_info->n_buffers = 2;
    d0_info->release = &ArrowAdapter::release_array;
    d0_info->buffers = new const void*[2];
    d0_info->buffers[0] = nullptr;
    d0_info->buffers[1] = malloc(sizeof(int64_t) * n);
    d0_info->n_children = 0;
    int64_t dom[] = {0, 1000, 1};
    std::memcpy((void*)d0_info->buffers[1], &dom, sizeof(int64_t) * n);

    return std::pair(
        std::move(arrow_schema),
        ArrowTable(std::move(col_info_array), std::move(col_info_schema)));
}

ArrowTable create_column_index_info() {
    // Create ArrowSchema for IndexColumnInfo
    auto col_info_schema = std::make_unique<ArrowSchema>();
    col_info_schema->format = "+s";
    col_info_schema->n_children = 1;
    col_info_schema->dictionary = nullptr;
    col_info_schema->release = &ArrowAdapter::release_schema;
    col_info_schema->children = new ArrowSchema*[col_info_schema->n_children];
    ArrowSchema* dim = col_info_schema->children[0] = new ArrowSchema;
    dim->format = "l";
    dim->name = "soma_dim_0";
    dim->n_children = 0;
    dim->dictionary = nullptr;
    dim->release = &ArrowAdapter::release_schema;

    // Create ArrowArray for IndexColumnInfo
    auto col_info_array = std::make_unique<ArrowArray>();
    col_info_array->length = 0;
    col_info_array->null_count = 0;
    col_info_array->offset = 0;
    col_info_array->n_buffers = 0;
    col_info_array->buffers = nullptr;
    col_info_array->n_children = 2;
    col_info_array->release = &ArrowAdapter::release_array;
    col_info_array->children = new ArrowArray*[1];
    int n = 3;
    auto d0_info = col_info_array->children[0] = new ArrowArray;
    d0_info->length = n;
    d0_info->null_count = 0;
    d0_info->offset = 0;
    d0_info->n_buffers = 2;
    d0_info->release = &ArrowAdapter::release_array;
    d0_info->buffers = new const void*[2];
    d0_info->buffers[0] = nullptr;
    d0_info->buffers[1] = malloc(sizeof(int64_t) * n);
    d0_info->n_children = 0;
    int64_t dom[] = {0, 1000, 1, 0, 2147483646};
    std::memcpy((void*)d0_info->buffers[1], &dom, sizeof(int64_t) * n);

    return ArrowTable(std::move(col_info_array), std::move(col_info_schema));
}
}  // namespace helper
