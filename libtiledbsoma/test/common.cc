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

std::pair<std::shared_ptr<ArrowSchema>, ColumnIndexInfo> create_arrow_schema() {
    // Create ArrowSchema
    auto arrow_schema = std::make_shared<ArrowSchema>();
    arrow_schema->format = "+s";
    arrow_schema->n_children = 2;
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->children = new ArrowSchema*[arrow_schema->n_children];

    ArrowSchema* dim = nullptr;
    dim = arrow_schema->children[0] = new ArrowSchema;
    dim->format = "l";
    dim->name = "d0";
    dim->n_children = 0;
    dim->dictionary = nullptr;
    dim->release = &ArrowAdapter::release_schema;

    ArrowSchema* attr = nullptr;
    attr = arrow_schema->children[1] = new ArrowSchema;
    attr->format = "l";
    attr->name = "a0";
    attr->n_children = 0;
    attr->dictionary = nullptr;
    attr->release = &ArrowAdapter::release_schema;

    // Create array for index columns
    std::vector<std::string> index_column_names = {"d0"};

    auto domains = std::make_shared<ArrowArray>();
    domains->length = 0;
    domains->null_count = 0;
    domains->offset = 0;
    domains->n_buffers = 0;
    domains->buffers = nullptr;
    domains->n_children = 2;
    domains->release = &ArrowAdapter::release_array;
    domains->children = new ArrowArray*[1];

    auto d0_domain = domains->children[0] = new ArrowArray;
    d0_domain->length = 2;
    d0_domain->null_count = 0;
    d0_domain->offset = 0;
    d0_domain->n_buffers = 2;
    d0_domain->release = &ArrowAdapter::release_array;
    d0_domain->buffers = new const void*[2];
    d0_domain->buffers[0] = nullptr;
    d0_domain->buffers[1] = malloc(sizeof(int64_t) * 2);
    d0_domain->n_children = 0;
    int64_t dom[] = {0, 1000};
    std::memcpy((void*)d0_domain->buffers[1], &dom, sizeof(int64_t) * 2);

    auto tiles = std::make_shared<ArrowArray>();
    tiles->length = 0;
    tiles->null_count = 0;
    tiles->offset = 0;
    tiles->n_buffers = 0;
    tiles->buffers = nullptr;
    tiles->n_children = 2;
    tiles->release = &ArrowAdapter::release_array;
    tiles->children = new ArrowArray*[1];

    ArrowArray* d0_tile = tiles->children[0] = new ArrowArray;
    d0_tile->length = 1;
    d0_tile->null_count = 0;
    d0_tile->offset = 0;
    d0_tile->n_buffers = 2;
    d0_tile->release = &ArrowAdapter::release_array;
    d0_tile->buffers = new const void*[2];
    d0_tile->buffers[0] = nullptr;
    d0_tile->buffers[1] = malloc(sizeof(int64_t));
    d0_tile->n_children = 0;
    int64_t tile = 1;
    std::memcpy((void*)d0_tile->buffers[1], &tile, sizeof(int64_t));

    ColumnIndexInfo index_columns_info = std::tuple(
        index_column_names, domains, tiles);

    return std::pair(arrow_schema, index_columns_info);
}
}  // namespace helper