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
#include "utils/logger.h"

namespace helper {

// This non-obvious number is:
// * Something that fits into signed 32-bit integer for R-friendliness;
// * Is a comfortable tile-extent distance away from 2^31-1 for default
//   core tile extent. (Using 2^31-1 exactly would result in a core
//   array-creation error.)
const int CORE_DOMAIN_MAX = 2147483646;

static std::unique_ptr<ArrowSchema> _create_index_cols_info_schema(
    const std::vector<DimInfo>& dim_infos);
static std::unique_ptr<ArrowArray> _create_index_cols_info_array(
    const std::vector<DimInfo>& dim_infos);

static std::string tdb_to_arrow_type(tiledb_datatype_t tiledb_dtype);

// Notes:
//
// * This is multi-purpose code used for generic SOMASparseNDArray,
//   SOMADenseNDArray, and SOMADataFrame
// * The NDArrays always have:
//   o int64 dims soma_dim_0, soma_dim_1, ..., soma_dim_n
//   o a single attr of varying numeric type -- float32, float32, int16, etc.
// * Default-indexed (i.e. what almost everyone uses) SOMADataFrame has
//   o int64 dims soma_joinid
//   o arbitrary number of attrs of arbitrary type -- int, float, bool, string,
//     categorical-of-string, you name it
// * But SOMADataFrame, in the SOMA spec, just needs to have soma_joinid
//   present, as a dim or attr.
//   o soma_joinid can be a dim, along with others
//   o soma_joinid can be an attr only, with something else being the dim
//
// Given this context, create_arrow_schema_and_index_columns is a factory
// for creating schema information for ND arrays and default-indexed
// dataframes.
//
// * This returns a pair of ArrowSchema and ArrowTable
// * ArrowTable, in turn, is a pair of ArrowArray and ArrowSchema
// * So this is
//    o ArrowSchema -- schema for the array
//      - ArrowArray  -- information about the tile extent, domain, and
//        current domain for the index columns
//      - ArrowSchema -- information about the datatypes for the index columns
//    o Confusingly, there's even more nesting: for an n-attr dataframe, there
//      is an ArrowSchema with its n_children being n -- and each of those
//      n entries in its children[] array are also of type ArrowSchema
// * The data structures conform to arrow_adapter
// * The Python and R bindings prepare similar Arrow information when
//   passing a create-array request to libtiledbsoma.

// Create ArrowSchema for the entire SOMAArray -- dims and attrs both -- as well
// as index-column info
std::pair<std::unique_ptr<ArrowSchema>, ArrowTable>
create_arrow_schema_and_index_columns(
    const std::vector<DimInfo>& dim_infos,
    const std::vector<AttrInfo>& attr_infos) {
    int ndim = dim_infos.size();
    int nattr = attr_infos.size();

    auto arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = "+s";
    arrow_schema->n_children = ndim + nattr;  // non-leaf node
    arrow_schema->dictionary = nullptr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->children = new ArrowSchema*[arrow_schema->n_children];

    for (int i = 0; i < ndim; i++) {
        const DimInfo& info = dim_infos[i];
        ArrowSchema* dim = new ArrowSchema;
        dim->name = strdup(info.name.c_str());
        dim->format = strdup(tdb_to_arrow_type(info.tiledb_datatype).c_str());
        dim->n_children = 0;  // leaf node
        dim->dictionary = nullptr;
        dim->release = &ArrowAdapter::release_schema;
        arrow_schema->children[i] = dim;
    }

    for (int i = 0; i < nattr; i++) {
        const AttrInfo& info = attr_infos[i];
        ArrowSchema* attr = new ArrowSchema;
        attr->name = strdup(info.name.c_str());
        attr->format = strdup(tdb_to_arrow_type(info.tiledb_datatype).c_str());
        attr->n_children = 0;  // leaf node
        attr->flags = 0;
        attr->dictionary = nullptr;
        attr->release = &ArrowAdapter::release_schema;
        arrow_schema->children[ndim + i] = attr;
    }

    auto index_cols_info_schema = _create_index_cols_info_schema(dim_infos);
    auto index_cols_info_array = _create_index_cols_info_array(dim_infos);

    return std::pair(
        std::move(arrow_schema),
        ArrowTable(
            std::move(index_cols_info_array),
            std::move(index_cols_info_schema)));
}

// Create index-column info only, no schema involving the attrs
ArrowTable create_column_index_info(const std::vector<DimInfo>& dim_infos) {
    for (auto info : dim_infos) {
        LOG_DEBUG(fmt::format(
            "create_column_index_info name={} type={} dim_max={} ucd={}",
            info.name,
            info.tiledb_datatype,
            info.dim_max,
            info.use_current_domain));
    }

    auto index_cols_info_schema = _create_index_cols_info_schema(dim_infos);
    auto index_cols_info_array = _create_index_cols_info_array(dim_infos);

    return ArrowTable(
        std::move(index_cols_info_array), std::move(index_cols_info_schema));
}

static std::unique_ptr<ArrowSchema> _create_index_cols_info_schema(
    const std::vector<DimInfo>& dim_infos) {
    int ndim = dim_infos.size();

    auto index_cols_info_schema = std::make_unique<ArrowSchema>();
    index_cols_info_schema->format = "+s";
    index_cols_info_schema->n_children = ndim;  // non-leaf node
    index_cols_info_schema->dictionary = nullptr;
    index_cols_info_schema->release = &ArrowAdapter::release_schema;
    index_cols_info_schema
        ->children = new ArrowSchema*[index_cols_info_schema->n_children];

    for (int i = 0; i < ndim; i++) {
        const DimInfo& info = dim_infos[i];
        ArrowSchema* dim_schema = new ArrowSchema;
        dim_schema->format = strdup(
            tdb_to_arrow_type(info.tiledb_datatype).c_str());
        dim_schema->name = strdup(info.name.c_str());
        dim_schema->n_children = 0;  // leaf node
        dim_schema->dictionary = nullptr;
        dim_schema->release = &ArrowAdapter::release_schema;
        index_cols_info_schema->children[i] = dim_schema;
    }

    return index_cols_info_schema;
}

static std::unique_ptr<ArrowArray> _create_index_cols_info_array(
    const std::vector<DimInfo>& dim_infos) {
    int ndim = dim_infos.size();

    auto index_cols_info_array = std::make_unique<ArrowArray>();
    index_cols_info_array->length = 0;
    index_cols_info_array->null_count = 0;
    index_cols_info_array->offset = 0;
    index_cols_info_array->n_buffers = 0;
    index_cols_info_array->buffers = nullptr;
    index_cols_info_array->n_children = ndim;
    index_cols_info_array->release = &ArrowAdapter::release_array;
    index_cols_info_array->children = new ArrowArray*[ndim];

    for (int i = 0; i < ndim; i++) {
        const DimInfo& info = dim_infos[i];

        int n = info.use_current_domain ? 5 : 3;

        auto dim_array = new ArrowArray;
        dim_array->length = n;
        dim_array->null_count = 0;
        dim_array->offset = 0;
        dim_array->n_buffers = 2;
        dim_array->release = &ArrowAdapter::release_array;
        dim_array->buffers = new const void*[2];
        dim_array->buffers[0] = nullptr;
        dim_array->n_children = 0;  // leaf node
        index_cols_info_array->children[i] = dim_array;

        // The full user-level SOMA API supports many more index types.
        // Here we support enough types to verify we've got variant-indexed
        // SOMADataFrame objects baseline-tested in C++, then defer exhaustive
        // loop-over-all-datatypes handling to Python and R.
        if (info.tiledb_datatype == TILEDB_INT64) {
            size_t nbytes = n * sizeof(int64_t);
            dim_array->buffers[1] = malloc(nbytes);
            if (info.use_current_domain) {
                // domain big; current_domain small
                int64_t dom[] = {0, CORE_DOMAIN_MAX, 1, 0, info.dim_max};
                void* vsrc = (void*)&dom[0];
                std::memcpy((void*)dim_array->buffers[1], vsrc, nbytes);
            } else {
                // domain small; current_domain feature not being used
                int64_t dom[] = {0, info.dim_max, 1};
                void* vsrc = (void*)&dom[0];
                std::memcpy((void*)dim_array->buffers[1], vsrc, nbytes);
            }

        } else if (info.tiledb_datatype == TILEDB_UINT32) {
            size_t nbytes = n * sizeof(uint32_t);
            dim_array->buffers[1] = malloc(nbytes);
            if (info.use_current_domain) {
                // domain big; current_domain small
                uint32_t dom[] = {
                    0, (uint32_t)CORE_DOMAIN_MAX, 1, 0, (uint32_t)info.dim_max};
                void* vsrc = (void*)&dom[0];
                std::memcpy((void*)dim_array->buffers[1], vsrc, nbytes);
            } else {
                // domain small; current_domain feature not being used
                int64_t dom[] = {0, info.dim_max, 1};
                void* vsrc = (void*)&dom[0];
                std::memcpy((void*)dim_array->buffers[1], vsrc, nbytes);
            }

        } else if (info.tiledb_datatype == TILEDB_STRING_ASCII) {
            // Domain specification for strings is not supported in core. See
            // arrow_adapter for more info. We rely on arrow_adapter to also
            // handle this case.
            dim_array->buffers[1] = nullptr;

        } else {
            throw TileDBSOMAError(
                "Unsupported datatype encountered in unit test. You can add a "
                "new type if you like!");
        }
    }

    return index_cols_info_array;
}

// Just a keystroke-saver
static std::string tdb_to_arrow_type(tiledb_datatype_t tiledb_datatype) {
    return std::string(ArrowAdapter::to_arrow_format(tiledb_datatype));
}

std::string to_arrow_format(tiledb_datatype_t tiledb_datatype) {
    return std::string(ArrowAdapter::to_arrow_format(tiledb_datatype));
}

}  // namespace helper
