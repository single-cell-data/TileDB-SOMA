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
        if (attr.nullable()) {
            child->flags |= ARROW_FLAG_NULLABLE;
        } else {
            child->flags &= ~ARROW_FLAG_NULLABLE;
        }
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
            if (enmr.ordered()) {
                child->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
            } else {
                child->flags &= ~ARROW_FLAG_DICTIONARY_ORDERED;
            }
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

FilterList ArrowAdapter::_create_filter_list(
    std::string filters, std::shared_ptr<Context> ctx) {
    return ArrowAdapter::_create_filter_list(json::parse(filters), ctx);
}

FilterList ArrowAdapter::_create_filter_list(
    json filters, std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    for (auto filter : filters) {
        ArrowAdapter::_append_to_filter_list(filter_list, filter, ctx);
    }

    return filter_list;
}

FilterList ArrowAdapter::_create_attr_filter_list(
    std::string name,
    PlatformConfig platform_config,
    std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    if (platform_config.attrs.empty()) {
        filter_list.add_filter(Filter(*ctx, TILEDB_FILTER_ZSTD));
    } else {
        json attr_options = json::parse(platform_config.attrs);
        if (attr_options.find(name) != attr_options.end() &&
            attr_options[name].find("filters") != attr_options[name].end()) {
            filter_list = ArrowAdapter::_create_filter_list(
                attr_options[name]["filters"], ctx);
        } else {
            filter_list.add_filter(Filter(*ctx, TILEDB_FILTER_ZSTD));
        }
    }

    return filter_list;
}

FilterList ArrowAdapter::_create_dim_filter_list(
    std::string name,
    PlatformConfig platform_config,
    std::string soma_type,
    std::shared_ptr<Context> ctx) {
    FilterList filter_list(*ctx);

    if (platform_config.dims.empty()) {
        filter_list.add_filter(
            ArrowAdapter::_get_zstd_default(platform_config, soma_type, ctx));
    } else {
        json dim_options = json::parse(platform_config.dims);
        if (dim_options.find(name) != dim_options.end() &&
            dim_options[name].find("filters") != dim_options[name].end()) {
            filter_list = ArrowAdapter::_create_filter_list(
                dim_options[name]["filters"], ctx);
        } else {
            filter_list.add_filter(ArrowAdapter::_get_zstd_default(
                platform_config, soma_type, ctx));
        }
    }

    return filter_list;
}

Filter ArrowAdapter::_get_zstd_default(
    PlatformConfig platform_config,
    std::string soma_type,
    std::shared_ptr<Context> ctx) {
    Filter zstd_filter(*ctx, TILEDB_FILTER_ZSTD);
    if (soma_type == "SOMADataFrame") {
        zstd_filter.set_option(
            TILEDB_COMPRESSION_LEVEL, platform_config.dataframe_dim_zstd_level);
    } else if (soma_type == "SOMASparseNDArray") {
        zstd_filter.set_option(
            TILEDB_COMPRESSION_LEVEL,
            platform_config.sparse_nd_array_dim_zstd_level);
    } else if (soma_type == "SOMADenseNDArray") {
        zstd_filter.set_option(
            TILEDB_COMPRESSION_LEVEL,
            platform_config.dense_nd_array_dim_zstd_level);
    }
    return zstd_filter;
}

void ArrowAdapter::_append_to_filter_list(
    FilterList filter_list, json value, std::shared_ptr<Context> ctx) {
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
        {"NONE", TILEDB_FILTER_NONE},
    };

    try {
        if (value.is_string()) {
            filter_list.add_filter(Filter(*ctx, convert_filter.at(value)));
        } else {
            Filter filter(*ctx, convert_filter.at(value["name"]));
            for (auto& [key, value] : value.items()) {
                ArrowAdapter::_set_filter_option(filter, key, value);
            }
            filter_list.add_filter(filter);
        }
    } catch (std::out_of_range& e) {
        throw TileDBSOMAError(fmt::format(
            "Invalid filter {} passed to PlatformConfig", std::string(value)));
    }
}

void ArrowAdapter::_set_filter_option(
    Filter filter, std::string option_name, json value) {
    if (option_name == "name") {
        return;
    }

    std::map<std::string, tiledb_filter_option_t> convert_option = {
        {"COMPRESSION_LEVEL", TILEDB_COMPRESSION_LEVEL},
        {"BIT_WIDTH_MAX_WINDOW", TILEDB_BIT_WIDTH_MAX_WINDOW},
        {"POSITIVE_DELTA_MAX_WINDOW", TILEDB_POSITIVE_DELTA_MAX_WINDOW},
        {"SCALE_FLOAT_BYTEWIDTH", TILEDB_SCALE_FLOAT_BYTEWIDTH},
        {"SCALE_FLOAT_FACTOR", TILEDB_SCALE_FLOAT_FACTOR},
        {"SCALE_FLOAT_OFFSET", TILEDB_SCALE_FLOAT_OFFSET},
        {"WEBP_INPUT_FORMAT", TILEDB_WEBP_INPUT_FORMAT},
        {"WEBP_QUALITY", TILEDB_WEBP_QUALITY},
        {"WEBP_LOSSLESS", TILEDB_WEBP_LOSSLESS},
        {"COMPRESSION_REINTERPRET_DATATYPE",
         TILEDB_COMPRESSION_REINTERPRET_DATATYPE},
    };

    auto option = convert_option[option_name];
    switch (option) {
        case TILEDB_COMPRESSION_LEVEL:
            filter.set_option(option, value.get<int32_t>());
            break;
        case TILEDB_BIT_WIDTH_MAX_WINDOW:
        case TILEDB_POSITIVE_DELTA_MAX_WINDOW:
            filter.set_option(option, value.get<uint32_t>());
            break;
        case TILEDB_SCALE_FLOAT_BYTEWIDTH:
            filter.set_option(option, value.get<uint64_t>());
            break;
        case TILEDB_SCALE_FLOAT_FACTOR:
        case TILEDB_SCALE_FLOAT_OFFSET:
            filter.set_option(option, value.get<double>());
            break;
        case TILEDB_WEBP_QUALITY:
            filter.set_option(option, value.get<float>());
            break;
        case TILEDB_WEBP_INPUT_FORMAT:
        case TILEDB_WEBP_LOSSLESS:
        case TILEDB_COMPRESSION_REINTERPRET_DATATYPE:
            filter.set_option(option, value.get<uint8_t>());
            break;
        default:
            throw TileDBSOMAError(
                fmt::format("Invalid option {} passed to filter", option_name));
    }
}

Dimension ArrowAdapter::_create_dim(
    tiledb_datatype_t type,
    std::string name,
    const void* buff,
    std::shared_ptr<Context> ctx) {
    switch (type) {
        case TILEDB_STRING_ASCII:
            return Dimension::create(*ctx, name, type, nullptr, nullptr);
        case TILEDB_TIME_SEC:
        case TILEDB_TIME_MS:
        case TILEDB_TIME_US:
        case TILEDB_TIME_NS:
        case TILEDB_DATETIME_SEC:
        case TILEDB_DATETIME_MS:
        case TILEDB_DATETIME_US:
        case TILEDB_DATETIME_NS:
            return Dimension::create(
                *ctx, name, type, (uint64_t*)buff, (uint64_t*)buff + 2);
        case TILEDB_INT8:
            return ArrowAdapter::_create_dim_aux(ctx, name, (int8_t*)buff);
        case TILEDB_UINT8:
            return ArrowAdapter::_create_dim_aux(ctx, name, (uint8_t*)buff);
        case TILEDB_INT16:
            return ArrowAdapter::_create_dim_aux(ctx, name, (int16_t*)buff);
        case TILEDB_UINT16:
            return ArrowAdapter::_create_dim_aux(ctx, name, (uint16_t*)buff);
        case TILEDB_INT32:
            return ArrowAdapter::_create_dim_aux(ctx, name, (int32_t*)buff);
        case TILEDB_UINT32:
            return ArrowAdapter::_create_dim_aux(ctx, name, (uint32_t*)buff);
        case TILEDB_INT64:
            return ArrowAdapter::_create_dim_aux(ctx, name, (int64_t*)buff);
        case TILEDB_UINT64:
            return ArrowAdapter::_create_dim_aux(ctx, name, (uint64_t*)buff);
        case TILEDB_FLOAT32:
            return ArrowAdapter::_create_dim_aux(ctx, name, (float*)buff);
        case TILEDB_FLOAT64:
            return ArrowAdapter::_create_dim_aux(ctx, name, (double*)buff);
        default:
            throw TileDBSOMAError(fmt::format(
                "ArrowAdapter: Unsupported TileDB dimension: {} ",
                tiledb::impl::type_to_str(type)));
    }
}

tiledb_layout_t ArrowAdapter::_get_order(std::string order) {
    std::transform(
        order.begin(), order.end(), order.begin(), [](unsigned char c) {
            return std::tolower(c);
        });

    std::map<std::string, tiledb_layout_t> convert_order = {
        {"row-major", TILEDB_ROW_MAJOR},
        {"row", TILEDB_ROW_MAJOR},
        {"col-major", TILEDB_COL_MAJOR},
        {"column-major", TILEDB_COL_MAJOR},
        {"col", TILEDB_COL_MAJOR},
        {"hilbert", TILEDB_HILBERT},
        {"unordered", TILEDB_UNORDERED},
    };

    try {
        return convert_order[order];
    } catch (const std::out_of_range& e) {
        throw TileDBSOMAError(
            fmt::format("Invalid order {} passed to PlatformConfig", order));
    }
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

    schema.set_capacity(platform_config.capacity);

    if (!platform_config.offsets_filters.empty()) {
        schema.set_offsets_filter_list(ArrowAdapter::_create_filter_list(
            platform_config.offsets_filters, ctx));
    }

    if (!platform_config.validity_filters.empty()) {
        schema.set_validity_filter_list(ArrowAdapter::_create_filter_list(
            platform_config.validity_filters, ctx));
    }

    schema.set_allows_dups(platform_config.allows_duplicates);

    if (platform_config.tile_order) {
        schema.set_tile_order(
            ArrowAdapter::_get_order(*platform_config.tile_order));
    }

    if (platform_config.cell_order) {
        schema.set_cell_order(
            ArrowAdapter::_get_order(*platform_config.cell_order));
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

                FilterList filter_list = ArrowAdapter::_create_dim_filter_list(
                    child->name, platform_config, soma_type, ctx);

                const void* buff = index_column_array->children[i]->buffers[1];
                auto dim = ArrowAdapter::_create_dim(
                    type, child->name, buff, ctx);
                dim.set_filter_list(filter_list);
                dims.insert({child->name, dim});
                isattr = false;
                break;
            }
        }

        if (isattr) {
            Attribute attr(*ctx, child->name, type);

            FilterList filter_list = ArrowAdapter::_create_attr_filter_list(
                child->name, platform_config, ctx);
            attr.set_filter_list(filter_list);

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
                LOG_DEBUG(fmt::format(
                    "[ArrowAdapter] dictionary for {} as {} {}",
                    child->name,
                    enmr_type,
                    enmr_format));
            }

            LOG_DEBUG(
                fmt::format("[ArrowAdapter] adding attribute {}", child->name));
            schema.add_attribute(attr);
        }
    }

    for (int64_t i = 0; i < index_column_schema->n_children; ++i) {
        LOG_DEBUG(fmt::format("[ArrowAdapter] child {}", i));
        auto col_name = index_column_schema->children[i]->name;
        domain.add_dimension(dims.at(col_name));
    }
    LOG_DEBUG(fmt::format("[ArrowAdapter] set_domain"));
    schema.set_domain(domain);

    LOG_DEBUG(fmt::format("[ArrowAdapter] check"));
    schema.check();

    LOG_DEBUG(fmt::format("[ArrowAdapter] returning"));
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
        for (size_t i = 0; i < column->size(); ++i) {
            array->null_count += column->validity()[i] == 0;
        }

        // Convert validity bytemap to a bitmap in place
        column->validity_to_bitmap();
        array->buffers[0] = column->validity().data();
    } else {
        schema->flags &= ~ARROW_FLAG_NULLABLE;  // as ArrowSchemaInitFromType
                                                // leads to NULLABLE set
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
