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

static std::unique_ptr<ArrowArray> _create_index_cols_info_array(
    const std::vector<DimInfo>& dim_infos);

// Core PR: https://github.com/TileDB-Inc/TileDB/pull/5303
bool have_dense_current_domain_support() {
    auto vers = tiledbsoma::version::embedded_version_triple();
    return std::get<0>(vers) >= 2 && std::get<1>(vers) >= 27;
}

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
    int nfield = ndim + nattr;

    std::vector<std::string> names(nfield);
    std::vector<tiledb_datatype_t> tiledb_datatypes(nfield);
    for (int i = 0; i < (int)ndim; i++) {
        const DimInfo& dim_info = dim_infos[i];
        names[i] = dim_info.name;
        tiledb_datatypes[i] = dim_info.tiledb_datatype;
    }
    for (int i = 0; i < (int)nattr; i++) {
        const AttrInfo& attr_info = attr_infos[i];
        names[ndim + i] = attr_info.name;
        tiledb_datatypes[ndim + i] = attr_info.tiledb_datatype;
    }
    auto arrow_schema = ArrowAdapter::make_arrow_schema(
        names, tiledb_datatypes);

    auto index_cols_info_schema = create_index_cols_info_schema(dim_infos);
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
            "create_column_index_info name={} type={} dim_max={}",
            info.name,
            tiledb::impl::to_str(info.tiledb_datatype),
            info.dim_max));
    }

    auto index_cols_info_schema = create_index_cols_info_schema(dim_infos);
    auto index_cols_info_array = _create_index_cols_info_array(dim_infos);

    return ArrowTable(
        std::move(index_cols_info_array), std::move(index_cols_info_schema));
}

std::unique_ptr<ArrowSchema> create_index_cols_info_schema(
    const std::vector<DimInfo>& dim_infos) {
    auto ndim = dim_infos.size();

    std::vector<std::string> names(ndim);
    std::vector<tiledb_datatype_t> tiledb_datatypes(ndim);

    for (int i = 0; i < static_cast<int>(ndim); i++) {
        const DimInfo& dim_info = dim_infos[i];
        names[i] = dim_info.name;
        tiledb_datatypes[i] = dim_info.tiledb_datatype;
    }

    auto schema = ArrowAdapter::make_arrow_schema(names, tiledb_datatypes);

    for (size_t i = 0; i < static_cast<size_t>(schema->n_children); ++i) {
        if (strcmp(schema->children[i]->name, "soma_geometry")) {
            nanoarrow::UniqueBuffer buffer;
            ArrowMetadataBuilderInit(buffer.get(), nullptr);
            ArrowMetadataBuilderAppend(
                buffer.get(),
                ArrowCharView("dtype"),
                ArrowCharView(
                    dim_infos[i].tiledb_datatype == TILEDB_GEOM_WKB ? "WKB" :
                                                                      "WKT"));
            ArrowSchemaSetMetadata(
                schema->children[i],
                std::string((char*)buffer->data, buffer->size_bytes).c_str());
        }
    }

    return schema;
}

static std::unique_ptr<ArrowArray> _create_index_cols_info_array(
    const std::vector<DimInfo>& dim_infos) {
    int ndim = dim_infos.size();

    auto index_cols_info_array = ArrowAdapter::make_arrow_array_parent(ndim);

    for (int i = 0; i < ndim; i++) {
        const DimInfo& info = dim_infos[i];
        ArrowArray* dim_array = nullptr;

        // The full user-level SOMA API supports many more index types.
        // Here we support enough types to verify we've got variant-indexed
        // SOMADataFrame objects baseline-tested in C++, then defer exhaustive
        // loop-over-all-datatypes handling to Python and R.
        if (info.tiledb_datatype == TILEDB_INT64) {
            // domain big; current_domain small
            std::vector<int64_t> dom({0, CORE_DOMAIN_MAX, 1, 0, info.dim_max});
            dim_array = ArrowAdapter::make_arrow_array_child(dom);
        } else if (info.tiledb_datatype == TILEDB_UINT32) {
            // domain big; current_domain small
            std::vector<uint32_t> dom(
                {0, (uint32_t)CORE_DOMAIN_MAX, 1, 0, (uint32_t)info.dim_max});
            dim_array = ArrowAdapter::make_arrow_array_child(dom);

        } else if (info.tiledb_datatype == TILEDB_STRING_ASCII) {
            // Domain specification for strings is not supported in core. See
            // arrow_adapter for more info. We rely on arrow_adapter to also
            // handle this case.
            std::vector<std::string> dom(
                {"", "", "", info.string_lo, info.string_hi});
            dim_array = ArrowAdapter::make_arrow_array_child_string(dom);
        } else if (info.tiledb_datatype == TILEDB_GEOM_WKB) {
            // No domain can be set for WKB. The domain will be set to the
            // individual spatial axes.
            dim_array = ArrowAdapter::make_arrow_array_child_binary();
        } else if (info.tiledb_datatype == TILEDB_FLOAT64) {
            // domain big; current_domain small
            std::vector<double_t> dom(
                {0, (double_t)CORE_DOMAIN_MAX, 1, 0, (double_t)info.dim_max});
            dim_array = ArrowAdapter::make_arrow_array_child(dom);
        }

        if (dim_array == nullptr) {
            throw TileDBSOMAError(
                "Unsupported datatype encountered in unit test. You can add a "
                "new type if you like!");
        }

        index_cols_info_array->children[i] = dim_array;
    }

    return index_cols_info_array;
}

}  // namespace helper
