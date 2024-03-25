/**
 * @file   arrow_adapter.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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
    schema->release = nullptr;

    for (int i = 0; i < schema->n_children; ++i) {
        struct ArrowSchema* child = schema->children[i];
        if (schema->name != nullptr) {
            free((void*)schema->name);
            schema->name = nullptr;
        }
        if (child->release != NULL) {
            child->release(child);
        }
        free(child);
    }
    free(schema->children);

    struct ArrowSchema* dict = schema->dictionary;
    if (dict != nullptr) {
        if (dict->format != nullptr) {
            free((void*)dict->format);
            dict->format = nullptr;
        }
        if (dict->release != nullptr) {
            delete dict;
            dict = nullptr;
        }
    }

    LOG_TRACE("[ArrowAdapter] release_schema");
}

void ArrowAdapter::release_array(struct ArrowArray* array) {
    auto arrow_buffer = static_cast<ArrowBuffer*>(array->private_data);

    LOG_TRACE(fmt::format(
        "[ArrowAdapter] release_array {} use_count={}",
        arrow_buffer->buffer_->name(),
        arrow_buffer->buffer_.use_count()));

    // Delete the ArrowBuffer, which was allocated with new.
    // If the ArrowBuffer.buffer_ shared_ptr is the last reference to the
    // underlying ColumnBuffer, the ColumnBuffer will be deleted.
    delete arrow_buffer;

    if (array->buffers != nullptr) {
        free(array->buffers);
    }

    struct ArrowArray* dict = array->dictionary;
    if (dict != nullptr) {
        if (dict->buffers != nullptr) {
            free(dict->buffers);
            dict->buffers = nullptr;
        }
        if (dict->release != nullptr) {
            delete dict;
            dict = nullptr;
        }
    }

    array->release = nullptr;
}

std::shared_ptr<ArrowSchema> ArrowAdapter::arrow_schema_from_tiledb_array(
    std::shared_ptr<Context> ctx, std::shared_ptr<Array> tiledb_array) {
    auto tiledb_schema = tiledb_array->schema();
    auto ndim = tiledb_schema.domain().ndim();
    auto nattr = tiledb_schema.attribute_num();

    std::shared_ptr<ArrowSchema> arrow_schema = std::make_shared<ArrowSchema>();
    arrow_schema->format = "+s";
    arrow_schema->n_children = ndim + nattr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->children = new ArrowSchema*[arrow_schema->n_children];

    ArrowSchema* child = nullptr;

    for (uint32_t i = 0; i < ndim; ++i) {
        auto dim = tiledb_schema.domain().dimension(i);
        child = arrow_schema->children[i] = new ArrowSchema;
        child->format = ArrowAdapter::to_arrow_format(dim.type()).data();
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
        child = arrow_schema->children[ndim + i] = new ArrowSchema;
        child->format = ArrowAdapter::to_arrow_format(attr.type()).data();
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
            auto dict = new ArrowSchema;
            dict->format = strdup(
                ArrowAdapter::to_arrow_format(enmr.type(), false).data());
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
            dst = new const void*[sz];
            std::memcpy((void*)dst, &src, sz);

            return std::pair(dst, data.size());
        }
        case TILEDB_INT8: {
            auto data = enmr.as_vector<int8_t>();
            return std::pair(
                ArrowAdapter::_fill_data_buffer(data, dst), data.size());
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

ArraySchema ArrowAdapter::tiledb_schema_from_arrow_schema(
    std::shared_ptr<Context> ctx,
    std::shared_ptr<ArrowSchema> arrow_schema,
    ColumnIndexInfo index_column_info,
    PlatformConfig platform_config) {
    auto [index_column_names, domains, extents] = index_column_info;

    std::cout << (platform_config["tiledb"]["create"]["allows_duplicates"] ?
                      "yes" :
                      "No")
              << std::endl;

    ArraySchema schema(*ctx, TILEDB_SPARSE);
    Domain domain(*ctx);

    std::map<std::string, Dimension> dims;

    for (int64_t sch_idx = 0; sch_idx < arrow_schema->n_children; ++sch_idx) {
        auto child = arrow_schema->children[sch_idx];
        auto type = ArrowAdapter::to_tiledb_format(child->format);

        auto idx_col_begin = index_column_names.begin();
        auto idx_col_end = index_column_names.end();
        auto idx_col_it = std::find(idx_col_begin, idx_col_end, child->name);

        if (idx_col_it != idx_col_end) {
            auto idx_col_idx = std::distance(idx_col_begin, idx_col_it);
            if (ArrowAdapter::_isvar(child->format)) {
                type = TILEDB_STRING_ASCII;
            }

            auto dim = Dimension::create(
                *ctx,
                child->name,
                type,
                type == TILEDB_STRING_ASCII ?
                    nullptr :
                    domains->children[idx_col_idx]->buffers[1],
                type == TILEDB_STRING_ASCII ?
                    nullptr :
                    extents->children[idx_col_idx]->buffers[1]);

            dims.insert({dim.name(), dim});
        } else {
            Attribute attr(*ctx, child->name, type);

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

    for (auto column_name : index_column_names)
        domain.add_dimension(dims.at(column_name));
    schema.set_domain(domain);

    schema.check();

    return schema;
}

ArrowTable ArrowAdapter::to_arrow(std::shared_ptr<ColumnBuffer> column) {
    std::shared_ptr<ArrowSchema> schema = std::make_shared<ArrowSchema>();
    std::shared_ptr<ArrowArray> array = std::make_shared<ArrowArray>();

    schema->format = to_arrow_format(column->type()).data();
    schema->name = column->name().data();
    schema->metadata = nullptr;
    schema->flags = 0;
    schema->n_children = 0;
    schema->children = nullptr;
    schema->dictionary = nullptr;
    schema->release = &release_schema;
    schema->private_data = nullptr;

    int n_buffers = column->is_var() ? 3 : 2;

    // Create an ArrowBuffer to manage the lifetime of `column`.
    // - `arrow_buffer` holds a shared_ptr to `column`, which
    // increments
    //   the use count and keeps the ColumnBuffer data alive.
    // - When the arrow array is released, `array->release()` is
    // called with
    //   `arrow_buffer` in `private_data`. `arrow_buffer` is
    //   deleted, which decrements the the `column` use count. When
    //   the `column` use count reaches 0, the ColumnBuffer data
    //   will be deleted.
    auto arrow_buffer = new ArrowBuffer(column);

    array->length = column->size();
    array->null_count = 0;
    array->offset = 0;
    array->n_buffers = n_buffers;
    array->n_children = 0;
    array->buffers = nullptr;
    array->children = nullptr;
    array->dictionary = nullptr;
    array->release = &release_array;
    array->private_data = (void*)arrow_buffer;

    LOG_TRACE(fmt::format(
        "[ArrowAdapter] create array name='{}' use_count={}",
        column->name(),
        column.use_count()));

    array->buffers = new const void*[n_buffers];
    assert(array->buffers != nullptr);
    array->buffers[0] = nullptr;                                   // validity
    array->buffers[n_buffers - 1] = column->data<void*>().data();  // data
    if (n_buffers == 3) {
        array->buffers[1] = column->offsets().data();  // offsets
    }

    if (column->is_nullable()) {
        schema->flags |= ARROW_FLAG_NULLABLE;

        // Count nulls
        for (auto v : column->validity()) {
            array->null_count += v == 0;
        }

        // Convert validity bytemap to a bitmap in place
        column->validity_to_bitmap();
        array->buffers[0] = column->validity().data();
    }
    if (column->is_ordered()) {
        schema->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
    }

    // Workaround to cast TILEDB_BOOL from uint8 to 1-bit Arrow boolean
    if (column->type() == TILEDB_BOOL) {
        column->data_to_bitmap();
    }

    if (column->has_enumeration()) {
        auto dict_sch = new ArrowSchema;
        auto dict_arr = new ArrowArray;

        auto enmr = column->get_enumeration_info();
        dict_sch->format = strdup(to_arrow_format(enmr->type(), false).data());
        dict_sch->name = nullptr;
        dict_sch->metadata = nullptr;
        dict_sch->flags = 0;
        dict_sch->n_children = 0;
        dict_sch->children = nullptr;
        dict_sch->dictionary = nullptr;
        dict_sch->release = &release_schema;
        dict_sch->private_data = nullptr;

        const int n_buf = ArrowAdapter::_isvar(dict_sch->format) ? 3 : 2;
        dict_arr->null_count = 0;
        dict_arr->offset = 0;
        dict_arr->n_buffers = n_buf;
        dict_arr->n_children = 0;
        dict_arr->buffers = nullptr;
        dict_arr->children = nullptr;
        dict_arr->dictionary = nullptr;
        dict_arr->release = &release_array;
        dict_arr->private_data = nullptr;

        dict_arr->buffers = new const void*[n_buf];
        dict_arr->buffers[0] = nullptr;  // validity: none here

        // TODO string types currently get the data and offset
        // buffers from ColumnBuffer::enum_offsets and
        // ColumnBuffer::enum_string which is retrieved via
        // ColumnBuffer::convert_enumeration. This may be refactored
        // to all use ColumnBuffer::get_enumeration_info. Note that
        // ColumnBuffer::has_enumeration may also be removed in a
        // future refactor as ColumnBuffer::get_enumeration_info
        // returns std::optional where std::nullopt indicates the
        // column does not contain enumerated values.
        if (enmr->type() == TILEDB_STRING_ASCII or
            enmr->type() == TILEDB_STRING_UTF8 or enmr->type() == TILEDB_CHAR) {
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

    return ArrowTable(array, schema);
}

bool ArrowAdapter::_isvar(const char* format) {
    if ((strcmp(format, "U") == 0) | (strcmp(format, "Z") == 0) |
        (strcmp(format, "u") == 0) | (strcmp(format, "z") == 0)) {
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

}  // namespace tiledbsoma