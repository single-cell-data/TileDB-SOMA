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
    LOG_DEBUG(fmt::format("[ArrowAdapter] release_schema for {}", schema->name));
    schema->release = nullptr;

    if (schema->name != nullptr) {
        free((void*)schema->name);
        schema->name = nullptr;
    }
    if (schema->format != nullptr) {
        free((void*)schema->format);
        schema->format = nullptr;
    }
    for (int i = 0; i < schema->n_children; ++i) {
        struct ArrowSchema* child = schema->children[i];
        if (child->name != nullptr) {
            free((void*)child->name);
            child->name = nullptr;
        }
        if (child->format != nullptr) {
            free((void*)child->format);
            child->format = nullptr;
        }
        if (child->release != NULL) {
            child->release(child);
        }
        free(child);
    }
    free(schema->children);

    struct ArrowSchema* dict = schema->dictionary;
    if (dict != nullptr) {
        if (dict->name != nullptr) {
            free((void*)dict->name);
            dict->name = nullptr;
        }
        if (dict->format != nullptr) {
            free((void*)dict->format);
            dict->format = nullptr;
        }
        if (dict->release != nullptr) {
            //delete dict;
            free(dict);
            dict = nullptr;
        }
    }
    LOG_TRACE("[ArrowAdapter] release_schema");
}

void ArrowAdapter::release_array(struct ArrowArray* array) {
    auto arrow_buffer = static_cast<ArrowBuffer*>(array->private_data);
    LOG_DEBUG(fmt::format("[ArrowAdapter] release_array for {} cnt {} var {} nullable {} enum {}",
                         arrow_buffer->buffer_->name(),
                         arrow_buffer->buffer_.use_count(),
                         arrow_buffer->buffer_->is_var(),
                         arrow_buffer->buffer_->is_nullable(),
                         arrow_buffer->buffer_->has_enumeration()
                         ));

    LOG_TRACE(fmt::format(
        "[ArrowAdapter] release_array {} use_count={}",
        arrow_buffer->buffer_->name(),
        arrow_buffer->buffer_.use_count()));

    // Delete the ArrowBuffer, which was allocated with new.
    // If the ArrowBuffer.buffer_ shared_ptr is the last reference to the
    // underlying ColumnBuffer, the ColumnBuffer will be deleted.
    delete arrow_buffer;

    if (array->buffers != nullptr) {
        delete[] array->buffers;
        array->buffers = nullptr;
    }

    if (array->n_children > 0) {
        for (int i = 0; i < array->n_children; ++i) {
            struct ArrowArray* child = array->children[i];
            if (child != nullptr) {
                release_array(child);
                free(child);
                child = nullptr;
            }
        }
        free(array->children);
        array->children = nullptr;
    }

    struct ArrowArray* dict = array->dictionary;
    if (dict != nullptr) {
        if (dict->buffers != nullptr) {
            free(dict->buffers);
            dict->buffers = nullptr;
        }
        if (dict->release != nullptr) {
            //delete dict;
            free(dict);
            dict = nullptr;
        }
    }
    array->release = nullptr;

}

std::unique_ptr<ArrowSchema> ArrowAdapter::arrow_schema_from_tiledb_array(
    std::shared_ptr<Context> ctx, std::shared_ptr<Array> tiledb_array) {
    auto tiledb_schema = tiledb_array->schema();
    auto ndim = tiledb_schema.domain().ndim();
    auto nattr = tiledb_schema.attribute_num();

    std::unique_ptr<ArrowSchema> arrow_schema = std::make_unique<ArrowSchema>();
    arrow_schema->format = "+s";
    arrow_schema->n_children = ndim + nattr;
    arrow_schema->release = &ArrowAdapter::release_schema;
    arrow_schema->children = (ArrowSchema**) malloc(arrow_schema->n_children * sizeof(ArrowSchema*)); //new ArrowSchema*[arrow_schema->n_children];

    ArrowSchema* child = nullptr;

    for (uint32_t i = 0; i < ndim; ++i) {
        auto dim = tiledb_schema.domain().dimension(i);
        child = arrow_schema->children[i] = (ArrowSchema*) malloc(sizeof(ArrowSchema)); //new ArrowSchema;
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
        child = arrow_schema->children[ndim + i] = (ArrowSchema*) malloc(sizeof(ArrowSchema)); //new ArrowSchema;
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
            if (enmr.type() == TILEDB_STRING_ASCII or
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
            dst = malloc(sz); //new const void*[sz];
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

bool ArrowAdapter::_isstr(const char* format) {
    if ((strcmp(format, "U") == 0) || (strcmp(format, "Z") == 0) ||
        (strcmp(format, "u") == 0) || (strcmp(format, "z") == 0)) {
        return true;
    }
    return false;
}

inline void exitIfError(const ArrowErrorCode ec, const std::string& msg) {
    if (ec != NANOARROW_OK)
        throw TileDBSOMAError(fmt::format(
                "ArrowAdapter: Arrow Error {} ", msg));
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
    exitIfError(ArrowSchemaSetName(sch, column->name().data()), "Bad schema name");
    exitIfError(ArrowSchemaAllocateChildren(sch, 0), "Bad schema children alloc");

#if 0
    schema->format = to_arrow_format(column->type()).data();
    schema->name = column->name().data();
    schema->metadata = nullptr;
    schema->flags = 0;
    schema->n_children = 0;
    schema->children = nullptr;
    schema->dictionary = nullptr;
#endif
    schema->release = &release_schema;
    schema->private_data = nullptr;

    int n_buffers = column->is_var() ? 3 : 2; // this will be 2 for enumerations and 3 for char vectors

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

    exitIfError(ArrowArrayInitFromType(arr, natype), "Bad array init");
    exitIfError(ArrowArrayAllocateChildren(arr, 0), "Bad array children alloc");
    array->length = column->size();

    LOG_DEBUG(fmt::format("[ArrowAdapter] column type {} name {} nbuf {} {} nullable {}",
                         to_arrow_format(column->type()).data(),
                         column->name().data(), n_buffers, array->n_buffers, column->is_nullable()));


#if 0
    array->null_count = 0;
    array->offset = 0;
    array->n_buffers = n_buffers;
    array->n_children = 0;
    array->buffers = nullptr;
    array->children = nullptr;
    array->dictionary = nullptr;
#endif
    array->release = &release_array;
    array->private_data = (void*)arrow_buffer;

    LOG_TRACE(fmt::format(
        "[ArrowAdapter] create array name='{}' use_count={}",
        column->name(),
        column.use_count()));

    array->buffers = (const void**) malloc(sizeof(void*) * n_buffers); //new const void*[n_buffers];
    assert(array->buffers != nullptr);
    array->buffers[0] = nullptr;                                   // validity addressed below
    array->buffers[n_buffers - 1] = column->data<void*>().data();  // data
    if (n_buffers == 3) {
        array->buffers[1] = column->offsets().data();  // offsets
    }

    if (column->is_nullable()) {
        schema->flags |= ARROW_FLAG_NULLABLE; // turns out it is also set by default

        // Count nulls
        for (auto v : column->validity()) {
            array->null_count += v == 0;
        }

        // Convert validity bytemap to a bitmap in place
        column->validity_to_bitmap();
        array->buffers[0] = column->validity().data();
    } else {
        schema->flags = 0;      // because ArrowSchemaInitFromType leads to NULLABLE set
    }

    if (column->is_ordered()) {
        schema->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
    }

    // Workaround to cast TILEDB_BOOL from uint8 to 1-bit Arrow boolean
    if (column->type() == TILEDB_BOOL) {
        column->data_to_bitmap();
    }

    if (column->has_enumeration()) {
        auto dict_sch = (ArrowSchema*) malloc(sizeof(ArrowSchema)); //new ArrowSchema;
        auto dict_arr = (ArrowArray*) malloc(sizeof(ArrowArray));  //new ArrowArray;

        auto enmr = column->get_enumeration_info();
        auto dcoltype = to_arrow_format(enmr->type(), false).data();
        auto dnatype = to_nanoarrow_type(dcoltype);

        exitIfError(ArrowSchemaInitFromType(dict_sch, dnatype), "Bad schema init");
        exitIfError(ArrowSchemaSetName(dict_sch, ""), "Bad schema name");
        exitIfError(ArrowSchemaAllocateChildren(dict_sch, 0), "Bad schema children alloc");
#if 0
        dict_sch->format = strdup(to_arrow_format(enmr->type(), false).data());
        dict_sch->name = nullptr;
        dict_sch->metadata = nullptr;
        dict_sch->flags = 0;
        dict_sch->n_children = 0;
        dict_sch->children = nullptr;
        dict_sch->dictionary = nullptr;
        dict_sch->release = &release_schema;
        dict_sch->private_data = nullptr;
#endif

        exitIfError(ArrowArrayInitFromType(dict_arr, dnatype), "Bad array init");
        exitIfError(ArrowArrayAllocateChildren(dict_arr, 0), "Bad array children alloc");
        const int n_buf = ArrowAdapter::_isstr(dict_sch->format) ? 3 : 2;
        dict_arr->buffers = (const void**) malloc(sizeof(void*) * n_buf); //new const void*[n_buf];
        dict_arr->buffers[0] = nullptr;  // validity: none here
        dict_arr->release = &release_array;
#if 0
        dict_arr->null_count = 0;
        dict_arr->offset = 0;
        dict_arr->n_buffers = n_buf;
        dict_arr->n_children = 0;
        dict_arr->buffers = nullptr;
        dict_arr->children = nullptr;
        dict_arr->dictionary = nullptr;
        dict_arr->private_data = nullptr;
#endif

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

    return std::pair(std::move(array), std::move(schema));
}

std::string_view ArrowAdapter::to_arrow_format(
    tiledb_datatype_t datatype, bool use_large) {
    switch (datatype) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
            return use_large ? "U" : "u";  // large because TileDB
                                           // uses 64bit offsets
        case TILEDB_CHAR:
        case TILEDB_BLOB:
            return use_large ? "Z" : "z";  // large because TileDB
                                           // uses 64bit offsets
        case TILEDB_BOOL:
            return "b";
        case TILEDB_INT32:
            return "i";
        case TILEDB_INT64:
            return "l";
        case TILEDB_FLOAT32:
            return "f";
        case TILEDB_FLOAT64:
            return "g";
        case TILEDB_INT8:
            return "c";
        case TILEDB_UINT8:
            return "C";
        case TILEDB_INT16:
            return "s";
        case TILEDB_UINT16:
            return "S";
        case TILEDB_UINT32:
            return "I";
        case TILEDB_UINT64:
            return "L";
        case TILEDB_TIME_SEC:
            return "tts";
        case TILEDB_TIME_MS:
            return "ttm";
        case TILEDB_TIME_US:
            return "ttu";
        case TILEDB_TIME_NS:
            return "ttn";
        case TILEDB_DATETIME_SEC:
            return "tss:";
        case TILEDB_DATETIME_MS:
            return "tsm:";
        case TILEDB_DATETIME_US:
            return "tsu:";
        case TILEDB_DATETIME_NS:
            return "tsn:";
        default:
            break;
    }
    throw TileDBSOMAError(fmt::format(
        "ArrowAdapter: Unsupported TileDB datatype: {} ",
        tiledb::impl::type_to_str(datatype)));
}

// FIXME: Add more types, maybe make it a map
enum ArrowType ArrowAdapter::to_nanoarrow_type(std::string_view sv) {
    if (sv == "i") 	 	 return NANOARROW_TYPE_INT32;
    else if (sv == "c")  return NANOARROW_TYPE_INT8;
    else if (sv == "C")  return NANOARROW_TYPE_UINT8;
    else if (sv == "s")  return NANOARROW_TYPE_INT16;
    else if (sv == "S")  return NANOARROW_TYPE_UINT16;
    else if (sv == "I")  return NANOARROW_TYPE_UINT32;
    else if (sv == "l")  return NANOARROW_TYPE_INT64;
    else if (sv == "L")  return NANOARROW_TYPE_UINT64;
    else if (sv == "f")  return NANOARROW_TYPE_FLOAT;
    else if (sv == "g")  return NANOARROW_TYPE_DOUBLE;
    else if (sv == "u")  return NANOARROW_TYPE_STRING;
    else if (sv == "U")  return NANOARROW_TYPE_LARGE_STRING;
    else if (sv == "b")  return NANOARROW_TYPE_BOOL;
    else throw TileDBSOMAError(fmt::format(
             "ArrowAdapter: Unsupported TileDB datatype string: {} ", sv));
}

}  // namespace tiledbsoma
