/**
 * @file   arrow_adapter.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022-2024 TileDB, Inc.
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
 * This file defines the ArrowAdapter class.
 */

#include "arrow_adapter.h"
#include "../soma/column_buffer.h"
#include "../utils/logger.h"

namespace tiledbsoma {

using namespace tiledb;

void ArrowAdapter::release_schema(struct ArrowSchema* schema) {
    if (schema->name != nullptr)
        LOG_DEBUG(
            fmt::format("[ArrowAdapter] release_schema for {}", schema->name));

    if (schema->name != nullptr) {
        LOG_TRACE("[ArrowAdapter] release_schema schema->name");
        free((void*)schema->name);
        schema->name = nullptr;
    }
    if (schema->format != nullptr) {
        LOG_TRACE("[ArrowAdapter] release_schema schema->format");
        free((void*)schema->format);
        schema->format = nullptr;
    }
    if (schema->metadata != nullptr) {
        LOG_TRACE("[ArrowAdapter] release_schema schema->metadata");
        free((void*)schema->metadata);
        schema->metadata = nullptr;
    }

    if (schema->children != nullptr) {
        for (auto i = 0; i < schema->n_children; i++) {
            if (schema->children[i] != nullptr) {
                if (schema->children[i]->release != nullptr) {
                    LOG_TRACE(fmt::format(
                        "[ArrowAdapter] release_schema schema->child {} "
                        "release",
                        i));
                    release_schema(schema->children[i]);
                }
                LOG_TRACE(fmt::format(
                    "[ArrowAdapter] release_schema schema->child {} free", i));
                free(schema->children[i]);
            }
        }
        LOG_TRACE("[ArrowAdapter] release_schema schema->children");
        free(schema->children);
        schema->children = nullptr;
    }

    if (schema->dictionary != nullptr) {
        if (schema->dictionary->release != nullptr) {
            LOG_TRACE("[ArrowAdapter] release_schema schema->dict release");
            release_schema(schema->dictionary);
        }
        LOG_TRACE("[ArrowAdapter] release_schema schema->dict free");
        free(schema->dictionary);
        schema->dictionary = nullptr;
    }

    schema->release = nullptr;
    LOG_TRACE("[ArrowAdapter] release_schema done");
}

void ArrowAdapter::release_array(struct ArrowArray* array) {
    auto arrow_buffer = static_cast<ArrowBuffer*>(array->private_data);
    if (arrow_buffer != nullptr) {
        LOG_TRACE(fmt::format(
            "[ArrowAdapter] release_array {} use_count={}",
            arrow_buffer->buffer_->name(),
            arrow_buffer->buffer_.use_count()));

        // Delete the ArrowBuffer, which was allocated with new.
        // If the ArrowBuffer.buffer_ shared_ptr is the last reference to the
        // underlying ColumnBuffer, the ColumnBuffer will be deleted.
        delete arrow_buffer;
    }

    if (array->buffers != nullptr) {
        free(array->buffers);
        array->buffers = nullptr;
    }

    if (array->children != nullptr) {
        for (auto i = 0; i < array->n_children; i++) {
            if (array->children[i] != nullptr) {
                if (array->children[i]->release != nullptr) {
                    LOG_TRACE(fmt::format(
                        "[ArrowAdapter] release_schema array->child {} release",
                        i));
                    release_array(array->children[i]);
                }
                LOG_TRACE(fmt::format(
                    "[ArrowAdapter] release_schema array->child {} free", i));
                free(array->children[i]);
            }
        }
        LOG_TRACE("[ArrowAdapter] release_array array->children");
        free(array->children);
        array->children = nullptr;
    }

    if (array->dictionary != nullptr) {
        // TODO: This can lead to segfault on some data sets, could be caused
        //       by how we fill arrow data structures.  This should pass.
        // if (array->dictionary->release != nullptr) {
        //    LOG_TRACE("[ArrowAdapter] release_array array->dict release");
        //    release_array(array->dictionary);
        //}
        LOG_TRACE("[ArrowAdapter] release_array array->dict free");
        free(array->dictionary);
        array->dictionary = nullptr;
    }

    array->release = nullptr;
    LOG_TRACE(fmt::format("[ArrowAdapter] release_array done"));
}

std::unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_from_tiledb_array(
    std::shared_ptr<Context> ctx, std::shared_ptr<Array> tiledb_array) {
    auto tiledb_schema = tiledb_array->schema();
    auto ndim = tiledb_schema.domain().ndim();
    auto nattr = tiledb_schema.attribute_num();

    std::unique_ptr<ArrowSchema> arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = strdup("+s");
    arrow_schema->n_children = ndim + nattr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->children = (ArrowSchema**)malloc(
        arrow_schema->n_children * sizeof(ArrowSchema*));

    ArrowSchema* child = nullptr;

    for (uint32_t i = 0; i < ndim; ++i) {
        auto dim = tiledb_schema.domain().dimension(i);
        child = arrow_schema->children[i] = (ArrowSchema*)malloc(
            sizeof(ArrowSchema));
        child->format = strdup(
            ArrowAdapter::to_arrow_format(dim.type()).data());
        child->name = strdup(dim.name().c_str());
        child->metadata = nullptr;
        child->flags = 0;
        child->n_children = 0;
        child->dictionary = nullptr;
        child->children = nullptr;
        child->release = &ArrowAdapter::release_schema;
    }

    for (uint32_t i = 0; i < nattr; ++i) {
        auto attr = tiledb_schema.attribute(i);
        child = arrow_schema->children[ndim + i] = (ArrowSchema*)malloc(
            sizeof(ArrowSchema));
        child->format = strdup(
            ArrowAdapter::to_arrow_format(attr.type()).data());
        child->name = strdup(attr.name().c_str());
        child->metadata = nullptr;
        child->flags = attr.nullable() ? ARROW_FLAG_NULLABLE : 0;
        child->n_children = 0;
        child->children = nullptr;
        child->dictionary = nullptr;

        auto enmr_name = AttributeExperimental::get_enumeration_name(
            *ctx, attr);
        if (enmr_name.has_value()) {
            auto enmr = ArrayExperimental::get_enumeration(
                *ctx, *tiledb_array, attr.name());
            auto dict = (ArrowSchema*)malloc(sizeof(ArrowSchema));
            dict->format = strdup(
                ArrowAdapter::to_arrow_format(enmr.type(), false).data());
            if (enmr.type() == TILEDB_STRING_ASCII ||
                enmr.type() == TILEDB_CHAR) {
                dict->format = strdup("z");
            } else {
                dict->format = strdup(
                    ArrowAdapter::to_arrow_format(enmr.type(), false).data());
            }
            dict->name = strdup(enmr.name().c_str());
            dict->metadata = nullptr;
            dict->flags = 0;
            dict->n_children = 0;
            dict->children = nullptr;
            dict->dictionary = nullptr;
            dict->release = &ArrowAdapter::release_schema;
            dict->private_data = nullptr;
            child->dictionary = dict;
        }
        child->release = &ArrowAdapter::release_schema;
    }

    return arrow_schema;
}

ArraySchema ArrowAdapter::tiledb_schema_from_arrow_schema(
    std::shared_ptr<Context> ctx,
    std::unique_ptr<ArrowSchema> arrow_schema,
    ArrowTable index_column_info,
    std::string soma_type,
    bool is_sparse,
    PlatformConfig platform_config) {
    auto index_column_array = std::move(index_column_info.first);
    auto index_column_schema = std::move(index_column_info.second);

    ArraySchema schema(*ctx, is_sparse ? TILEDB_SPARSE : TILEDB_DENSE);
    Domain domain(*ctx);

    std::map<std::string, tiledb_filter_type_t> convert_filter = {
        {"GZIP", TILEDB_FILTER_GZIP},
        {"ZSTD", TILEDB_FILTER_ZSTD},
        {"LZ4", TILEDB_FILTER_LZ4},
        {"BZIP2", TILEDB_FILTER_BZIP2},
        {"RLE", TILEDB_FILTER_RLE},
        {"DELTA", TILEDB_FILTER_DELTA},
        {"DOUBLE_DELTA", TILEDB_FILTER_DOUBLE_DELTA},
        {"BIT_WIDTH_REDUCTION", TILEDB_FILTER_BIT_WIDTH_REDUCTION},
        {"BITSHUFFLE", TILEDB_FILTER_BITSHUFFLE},
        {"BYTESHUFFLE", TILEDB_FILTER_BYTESHUFFLE},
        {"POSITIVE_DELTA", TILEDB_FILTER_POSITIVE_DELTA},
        {"CHECKSUM_MD5", TILEDB_FILTER_CHECKSUM_MD5},
        {"CHECKSUM_SHA256", TILEDB_FILTER_CHECKSUM_SHA256},
        {"DICTIONARY_ENCODING", TILEDB_FILTER_DICTIONARY},
        {"SCALE_FLOAT", TILEDB_FILTER_SCALE_FLOAT},
        {"XOR", TILEDB_FILTER_XOR},
        {"WEBP", TILEDB_FILTER_WEBP},
        {"NOOP", TILEDB_FILTER_NONE},
    };

    schema.set_capacity(platform_config.capacity);

    if (!platform_config.offsets_filters.empty()) {
        FilterList offset_filter_list(*ctx);
        json offsets_filters = json::parse(platform_config.offsets_filters);
        for (auto offset : offsets_filters) {
            try {
                offset_filter_list.add_filter(
                    Filter(*ctx, convert_filter.at(offset)));
            } catch (std::out_of_range& e) {
                throw TileDBSOMAError(fmt::format(
                    "Invalid offset filter {} passed to PlatformConfig",
                    offset));
            }
        }
        schema.set_offsets_filter_list(offset_filter_list);
    }

    if (!platform_config.validity_filters.empty()) {
        FilterList offset_filter_list(*ctx);
        json validity_filters = json::parse(platform_config.validity_filters);
        for (auto offset : validity_filters) {
            try {
                offset_filter_list.add_filter(
                    Filter(*ctx, convert_filter.at(offset)));
            } catch (std::out_of_range& e) {
                throw TileDBSOMAError(fmt::format(
                    "Invalid offset filter {} passed to PlatformConfig",
                    offset));
            }
        }
        schema.set_validity_filter_list(offset_filter_list);
    }

    schema.set_allows_dups(platform_config.allows_duplicates);

    if (platform_config.tile_order) {
        auto tile_order = *platform_config.tile_order;
        std::transform(
            tile_order.begin(),
            tile_order.end(),
            tile_order.begin(),
            [](unsigned char c) { return std::tolower(c); });

        if (tile_order == "row-order" or tile_order == "row") {
            schema.set_tile_order(TILEDB_ROW_MAJOR);
        } else if (tile_order == "col-order" or tile_order == "col") {
            schema.set_tile_order(TILEDB_COL_MAJOR);
        } else {
            throw TileDBSOMAError(fmt::format(
                "Invalid tile order {} passed to PlatformConfig", tile_order));
        }
    }

    if (platform_config.cell_order) {
        auto cell_order = *platform_config.cell_order;
        std::transform(
            cell_order.begin(),
            cell_order.end(),
            cell_order.begin(),
            [](unsigned char c) { return std::tolower(c); });

        if (cell_order == "row-major" or cell_order == "row") {
            schema.set_cell_order(TILEDB_ROW_MAJOR);
        } else if (
            cell_order == "col-major" or cell_order == "column-major" or
            cell_order == "col") {
            schema.set_cell_order(TILEDB_COL_MAJOR);
        } else if (cell_order == "hilbert") {
            schema.set_cell_order(TILEDB_HILBERT);
        } else {
            throw TileDBSOMAError(fmt::format(
                "Invalid cell order {} passed to PlatformConfig", cell_order));
        }
    }

    std::map<std::string, Dimension> dims;

    for (int64_t sch_idx = 0; sch_idx < arrow_schema->n_children; ++sch_idx) {
        auto child = arrow_schema->children[sch_idx];
        auto type = ArrowAdapter::to_tiledb_format(child->format);

        bool isattr = true;

        for (int64_t i = 0; i < index_column_schema->n_children; ++i) {
            auto col_name = index_column_schema->children[i]->name;
            if (strcmp(child->name, col_name) == 0) {
                if (ArrowAdapter::_isvar(child->format)) {
                    type = TILEDB_STRING_ASCII;
                }

                FilterList dim_filterlist(*ctx);
                Filter dim_zstd_filter(*ctx, TILEDB_FILTER_ZSTD);
                if (soma_type == "SOMADataFrame") {
                    dim_zstd_filter.set_option(
                        TILEDB_COMPRESSION_LEVEL,
                        platform_config.dataframe_dim_zstd_level);
                } else if (soma_type == "SOMASparseNDArray") {
                    dim_zstd_filter.set_option(
                        TILEDB_COMPRESSION_LEVEL,
                        platform_config.sparse_nd_array_dim_zstd_level);
                } else if (soma_type == "SOMADenseNDArray") {
                    dim_zstd_filter.set_option(
                        TILEDB_COMPRESSION_LEVEL,
                        platform_config.dense_nd_array_dim_zstd_level);
                }
                dim_filterlist.add_filter(dim_zstd_filter);

                const void* buff = index_column_array->children[i]->buffers[1];

                switch (type) {
                    case TILEDB_STRING_ASCII: {
                        auto dim = Dimension::create(
                            *ctx, child->name, type, nullptr, nullptr);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_TIME_SEC:
                    case TILEDB_TIME_MS:
                    case TILEDB_TIME_US:
                    case TILEDB_TIME_NS:
                    case TILEDB_DATETIME_SEC:
                    case TILEDB_DATETIME_MS:
                    case TILEDB_DATETIME_US:
                    case TILEDB_DATETIME_NS: {
                        auto datetime_buff = (uint64_t*)buff;
                        auto dim = Dimension::create(
                            *ctx,
                            child->name,
                            type,
                            datetime_buff,
                            datetime_buff + 2);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_INT8: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (int8_t*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_UINT8: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (uint8_t*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_INT16: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (int16_t*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_UINT16: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (uint16_t*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_INT32: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (int32_t*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_UINT32: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (uint32_t*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_INT64: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (int64_t*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_UINT64: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (uint64_t*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_FLOAT32: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (float*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    case TILEDB_FLOAT64: {
                        auto dim = ArrowAdapter::_create_dim(
                            *ctx, child->name, (double*)buff);
                        dim.set_filter_list(dim_filterlist);
                        dims.insert({child->name, dim});
                        break;
                    }
                    default:
                        throw TileDBSOMAError(fmt::format(
                            "ArrowAdapter: Unsupported TileDB dimension: {} ",
                            tiledb::impl::type_to_str(type)));
                }
                isattr = false;
                break;
            }
        }
        if (isattr) {
            Attribute attr(*ctx, child->name, type);

            FilterList attr_filterlist(*ctx);
            Filter attr_zstd_filter(*ctx, TILEDB_FILTER_ZSTD);
            attr_filterlist.add_filter(attr_zstd_filter);
            attr.set_filter_list(attr_filterlist);

            if (child->flags & ARROW_FLAG_NULLABLE) {
                attr.set_nullable(true);
            }

            if (ArrowAdapter::_isvar(child->format)) {
                attr.set_cell_val_num(TILEDB_VAR_NUM);
            }

            if (child->dictionary != nullptr) {
                auto enmr_format = child->dictionary->format;
                auto enmr_type = ArrowAdapter::to_tiledb_format(enmr_format);
                auto enmr = Enumeration::create_empty(
                    *ctx,
                    child->name,
                    enmr_type,
                    ArrowAdapter::_isvar(enmr_format) ? TILEDB_VAR_NUM : 1,
                    child->flags & ARROW_FLAG_DICTIONARY_ORDERED);
                ArraySchemaExperimental::add_enumeration(*ctx, schema, enmr);
                AttributeExperimental::set_enumeration_name(
                    *ctx, attr, child->name);
            }

            schema.add_attribute(attr);
        }
    }

    for (int64_t i = 0; i < index_column_schema->n_children; ++i) {
        auto col_name = index_column_schema->children[i]->name;
        domain.add_dimension(dims.at(col_name));
    }
    schema.set_domain(domain);

    schema.check();

    return schema;
}

std::pair<const void*, std::size_t> ArrowAdapter::_get_data_and_length(
    Enumeration& enmr, const void* dst) {
    switch (enmr.type()) {
        case TILEDB_BOOL: {
            // We must handle this specially because vector<bool> does
            // not store elements contiguously in memory
            auto data = enmr.as_vector<bool>();

            // Represent the Boolean vector with, at most, the last two
            // bits. In Arrow, Boolean values are LSB packed
            uint8_t src = 0;
            for (size_t i = 0; i < data.size(); ++i)
                src |= (data[i] << i);

            // Allocate a single byte to copy the bits into
            size_t sz = 1;
            dst = malloc(sz);
            std::memcpy((void*)dst, &src, sz);

            return std::pair(dst, data.size());
        }
        case TILEDB_INT8: {
            auto data = enmr.as_vector<int8_t>();
            return std::pair(_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_UINT8: {
            auto data = enmr.as_vector<uint8_t>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_INT16: {
            auto data = enmr.as_vector<int16_t>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_UINT16: {
            auto data = enmr.as_vector<uint16_t>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_INT32: {
            auto data = enmr.as_vector<int32_t>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_UINT32: {
            auto data = enmr.as_vector<uint32_t>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_INT64: {
            auto data = enmr.as_vector<int64_t>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_UINT64: {
            auto data = enmr.as_vector<uint64_t>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_FLOAT32: {
            auto data = enmr.as_vector<float>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        case TILEDB_FLOAT64: {
            auto data = enmr.as_vector<double>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
        }
        default:
            throw TileDBSOMAError(fmt::format(
                "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                tiledb::impl::type_to_str(enmr.type())));
    }
}

inline void exitIfError(const ArrowErrorCode ec, const std::string& msg) {
    if (ec != NANOARROW_OK)
        throw TileDBSOMAError(
            fmt::format("ArrowAdapter: Arrow Error {} ", msg));
}

std::pair<std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>>
ArrowAdapter::to_arrow(std::shared_ptr<ColumnBuffer> column) {
    std::unique_ptr<ArrowSchema> schema = std::make_unique<ArrowSchema>();
    std::unique_ptr<ArrowArray> array = std::make_unique<ArrowArray>();
    auto sch = schema.get();
    auto arr = array.get();

    auto coltype = to_arrow_format(column->type()).data();
    auto natype = to_nanoarrow_type(coltype);
    exitIfError(ArrowSchemaInitFromType(sch, natype), "Bad schema init");
    exitIfError(
        ArrowSchemaSetName(sch, column->name().data()), "Bad schema name");
    exitIfError(
        ArrowSchemaAllocateChildren(sch, 0), "Bad schema children alloc");
    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    schema->release = &release_schema;

    // this will be 3 for char vecs and 2 for enumerations
    int n_buffers = column->is_var() ? 3 : 2;

    // Create an ArrowBuffer to manage the lifetime of `column`.
    // - `arrow_buffer` holds shared_ptr to `column`, increments
    //   the use count and keeps the ColumnBuffer data alive.
    // - When the arrow array is released, `array->release()` is
    //   called with `arrow_buffer` in `private_data`.
    //   `arrow_buffer` is deleted, which decrements the the
    //   `column` use count. When the `column` use count reaches
    //   0, the ColumnBuffer data will be deleted.
    auto arrow_buffer = new ArrowBuffer(column);

    exitIfError(ArrowArrayInitFromType(arr, natype), "Bad array init");
    exitIfError(ArrowArrayAllocateChildren(arr, 0), "Bad array children alloc");
    array->length = column->size();

    LOG_TRACE(fmt::format(
        "[ArrowAdapter] column type {} name {} nbuf {} {} nullable {}",
        to_arrow_format(column->type()).data(),
        column->name().data(),
        n_buffers,
        array->n_buffers,
        column->is_nullable()));

    // After allocating and initializing via nanoarrow we
    // hook our custom release function in
    array->release = &release_array;
    if (array->private_data != nullptr) {  // as we use nanoarrow's init
        free(array->private_data);         // free what was allocated before
    }                                      // assigning our ArrowBuffer pointer
    array->private_data = (void*)arrow_buffer;

    LOG_TRACE(fmt::format(
        "[ArrowAdapter] create array name='{}' use_count={}",
        column->name(),
        column.use_count()));

    array->buffers = (const void**)malloc(sizeof(void*) * n_buffers);
    assert(array->buffers != nullptr);
    array->buffers[0] = nullptr;  // validity addressed below
    array->buffers[n_buffers - 1] = column->data<void*>().data();  // data
    if (n_buffers == 3) {
        array->buffers[1] = column->offsets().data();  // offsets
    }

    if (column->is_nullable()) {
        schema->flags |= ARROW_FLAG_NULLABLE;  // it is also set by default

        // Count nulls
        for (auto v : column->validity()) {
            array->null_count += v == 0;
        }

        // Convert validity bytemap to a bitmap in place
        column->validity_to_bitmap();
        array->buffers[0] = column->validity().data();
    } else {
        schema->flags = 0;  // as ArrowSchemaInitFromType leads to NULLABLE
                            // set
    }

    if (column->is_ordered()) {
        schema->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
    }

    // Workaround to cast TILEDB_BOOL from uint8 to 1-bit Arrow boolean
    if (column->type() == TILEDB_BOOL) {
        column->data_to_bitmap();
    }

    // Workaround for datetime
    if (column->type() == TILEDB_DATETIME_SEC ||
        column->type() == TILEDB_DATETIME_MS ||
        column->type() == TILEDB_DATETIME_NS) {
        free((void*)schema->format);  // free the 'storage' format
        schema->format = strdup(to_arrow_format(column->type()).data());
    }

    // Workaround for date
    if (column->type() == TILEDB_DATETIME_DAY) {
        free((void*)schema->format);  // free the 'storage' format
        schema->format = strdup(to_arrow_format(column->type()).data());
        // TODO: Put in ColumnBuffer
        size_t n = array->length;
        std::vector<int64_t> indata(n);
        std::memcpy(
            indata.data(), column->data<double>().data(), sizeof(int64_t) * n);
        std::vector<int32_t> vec(n);
        for (size_t i = 0; i < n; i++) {
            vec[i] = static_cast<int32_t>(indata[i]);
        }
        std::memcpy(
            (void*)array->buffers[n_buffers - 1],
            vec.data(),
            sizeof(int32_t) * n);
    }

    if (column->has_enumeration()) {
        auto dict_sch = (ArrowSchema*)malloc(sizeof(ArrowSchema));
        auto dict_arr = (ArrowArray*)malloc(sizeof(ArrowArray));

        auto enmr = column->get_enumeration_info();
        auto dcoltype = to_arrow_format(enmr->type(), false).data();
        auto dnatype = to_nanoarrow_type(dcoltype);

        exitIfError(
            ArrowSchemaInitFromType(dict_sch, dnatype), "Bad schema init");
        exitIfError(ArrowSchemaSetName(dict_sch, ""), "Bad schema name");
        exitIfError(
            ArrowSchemaAllocateChildren(dict_sch, 0),
            "Bad schema children alloc");
        dict_sch->release = &release_schema;

        exitIfError(
            ArrowArrayInitFromType(dict_arr, dnatype), "Bad array init");
        exitIfError(
            ArrowArrayAllocateChildren(dict_arr, 0),
            "Bad array children alloc");
        // hook up our custom release function
        dict_arr->release = &release_array;

        // TODO string types currently get the data and offset
        // buffers from ColumnBuffer::enum_offsets and
        // ColumnBuffer::enum_string which is retrieved via
        // ColumnBuffer::convert_enumeration. This may be refactored
        // to all use ColumnBuffer::get_enumeration_info. Note that
        // ColumnBuffer::has_enumeration may also be removed in a
        // future refactor as ColumnBuffer::get_enumeration_info
        // returns std::optional where std::nullopt indicates the
        // column does not contain enumerated values.
        if (enmr->type() == TILEDB_STRING_ASCII ||
            enmr->type() == TILEDB_STRING_UTF8 || enmr->type() == TILEDB_CHAR) {
            auto dict_vec = enmr->as_vector<std::string>();
            column->convert_enumeration();
            dict_arr->buffers[1] = column->enum_offsets().data();
            dict_arr->buffers[2] = column->enum_string().data();
            dict_arr->length = dict_vec.size();
        } else {
            auto [dict_data, dict_length] = _get_data_and_length(
                *enmr, dict_arr->buffers[1]);
            dict_arr->buffers[1] = dict_data;
            dict_arr->length = dict_length;
        }

        schema->dictionary = dict_sch;
        array->dictionary = dict_arr;
    }

    return std::pair(std::move(array), std::move(schema));
}

bool ArrowAdapter::_isvar(const char* format) {
    if ((strcmp(format, "U") == 0) || (strcmp(format, "Z") == 0) ||
        (strcmp(format, "u") == 0) || (strcmp(format, "z") == 0)) {
        return true;
    }
    return false;
}

std::string_view ArrowAdapter::to_arrow_format(
    tiledb_datatype_t tiledb_dtype, bool use_large) {
    auto u = use_large ? "U" : "u";
    auto z = use_large ? "Z" : "z";
    std::map<tiledb_datatype_t, std::string_view> _to_arrow_format_map = {
        {TILEDB_STRING_ASCII, u},     {TILEDB_CHAR, z},
        {TILEDB_STRING_UTF8, u},      {TILEDB_BLOB, z},
        {TILEDB_INT8, "c"},           {TILEDB_UINT8, "C"},
        {TILEDB_INT16, "s"},          {TILEDB_UINT16, "S"},
        {TILEDB_INT32, "i"},          {TILEDB_UINT32, "I"},
        {TILEDB_INT64, "l"},          {TILEDB_UINT64, "L"},
        {TILEDB_FLOAT32, "f"},        {TILEDB_FLOAT64, "g"},
        {TILEDB_BOOL, "b"},           {TILEDB_DATETIME_SEC, "tss:"},
        {TILEDB_DATETIME_MS, "tsm:"}, {TILEDB_DATETIME_US, "tsu:"},
        {TILEDB_DATETIME_NS, "tsn:"},
    };

    try {
        return _to_arrow_format_map.at(tiledb_dtype);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(fmt::format(
            "ArrowAdapter: Unsupported TileDB type: {} ",
            tiledb::impl::type_to_str(tiledb_dtype)));
    }
}

tiledb_datatype_t ArrowAdapter::to_tiledb_format(std::string_view arrow_dtype) {
    std::map<std::string_view, tiledb_datatype_t> _to_tiledb_format_map = {
        {"u", TILEDB_STRING_UTF8},    {"U", TILEDB_STRING_UTF8},
        {"z", TILEDB_CHAR},           {"Z", TILEDB_CHAR},
        {"c", TILEDB_INT8},           {"C", TILEDB_UINT8},
        {"s", TILEDB_INT16},          {"S", TILEDB_UINT16},
        {"i", TILEDB_INT32},          {"I", TILEDB_UINT32},
        {"l", TILEDB_INT64},          {"L", TILEDB_UINT64},
        {"f", TILEDB_FLOAT32},        {"g", TILEDB_FLOAT64},
        {"b", TILEDB_BOOL},           {"tss:", TILEDB_DATETIME_SEC},
        {"tsm:", TILEDB_DATETIME_MS}, {"tsu:", TILEDB_DATETIME_US},
        {"tsn:", TILEDB_DATETIME_NS},
    };

    try {
        return _to_tiledb_format_map.at(arrow_dtype);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(fmt::format(
            "ArrowAdapter: Unsupported Arrow type: {} ", arrow_dtype));
    }
}

// FIXME: Add more types, maybe make it a map
enum ArrowType ArrowAdapter::to_nanoarrow_type(std::string_view sv) {
    if (sv == "i")
        return NANOARROW_TYPE_INT32;
    else if (sv == "c")
        return NANOARROW_TYPE_INT8;
    else if (sv == "C")
        return NANOARROW_TYPE_UINT8;
    else if (sv == "s")
        return NANOARROW_TYPE_INT16;
    else if (sv == "S")
        return NANOARROW_TYPE_UINT16;
    else if (sv == "I")
        return NANOARROW_TYPE_UINT32;
    else if (sv == "l")
        return NANOARROW_TYPE_INT64;
    else if (sv == "L")
        return NANOARROW_TYPE_UINT64;
    else if (sv == "f")
        return NANOARROW_TYPE_FLOAT;
    else if (sv == "g")
        return NANOARROW_TYPE_DOUBLE;
    else if (sv == "u")
        return NANOARROW_TYPE_STRING;
    else if (sv == "U")
        return NANOARROW_TYPE_LARGE_STRING;
    else if (sv == "b")
        return NANOARROW_TYPE_BOOL;
    else if (sv == "tss:")
        return NANOARROW_TYPE_INT64;  // NB time resolution set indepedently
    else if (sv == "tsm:")
        return NANOARROW_TYPE_INT64;  // NB time resolution set indepedently
    else if (sv == "tsn:")
        return NANOARROW_TYPE_INT64;  // NB time resolution set indepedently
    else if (sv == "tsu:")
        return NANOARROW_TYPE_INT64;  // NB time resolution set indepedently
    else if (sv == "tdD")
        return NANOARROW_TYPE_INT32;  // R Date: fractional days since epoch
    else if (sv == "z")
        return NANOARROW_TYPE_BINARY;
    else if (sv == "Z")
        return NANOARROW_TYPE_LARGE_BINARY;
    else
        throw TileDBSOMAError(
            fmt::format("ArrowAdapter: Unsupported Arrow format: {} ", sv));
}

}  // namespace tiledbsoma
