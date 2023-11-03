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

std::pair<std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>>
ArrowAdapter::to_arrow(std::shared_ptr<ColumnBuffer> column, bool use_enum) {
    std::unique_ptr<ArrowSchema> schema = std::make_unique<ArrowSchema>();
    std::unique_ptr<ArrowArray> array = std::make_unique<ArrowArray>();

    schema->format = to_arrow_format(column->type()).data();  // mandatory
    schema->name = column->name().data();                     // optional
    schema->metadata = nullptr;                               // optional
    schema->flags = 0;                                        // optional
    schema->n_children = 0;                                   // mandatory
    schema->children = nullptr;                               // optional
    schema->dictionary = nullptr;                             // optional
    schema->release = &release_schema;                        // mandatory
    schema->private_data = nullptr;                           // optional

    int n_buffers = column->is_var() ? 3 : 2;

    // Create an ArrowBuffer to manage the lifetime of `column`.
    // - `arrow_buffer` holds a shared_ptr to `column`, which increments
    //   the use count and keeps the ColumnBuffer data alive.
    // - When the arrow array is released, `array->release()` is called with
    //   `arrow_buffer` in `private_data`. `arrow_buffer` is deleted, which
    //   decrements the the `column` use count. When the `column` use count
    //   reaches 0, the ColumnBuffer data will be deleted.
    auto arrow_buffer = new ArrowBuffer(column);

    array->length = column->size();             // mandatory
    array->null_count = 0;                      // mandatory
    array->offset = 0;                          // mandatory
    array->n_buffers = n_buffers;               // mandatory
    array->n_children = 0;                      // mandatory
    array->buffers = nullptr;                   // mandatory
    array->children = nullptr;                  // optional
    array->dictionary = nullptr;                // optional
    array->release = &release_array;            // mandatory
    array->private_data = (void*)arrow_buffer;  // mandatory

    LOG_TRACE(fmt::format(
        "[ArrowAdapter] create array name='{}' use_count={}",
        column->name(),
        column.use_count()));

    array->buffers = (const void**)malloc(sizeof(void*) * n_buffers);
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

    /* Workaround to cast TILEDB_BOOL from uint8 to 1-bit Arrow boolean. */
    if (column->type() == TILEDB_BOOL) {
        column->data_to_bitmap();
    }

    // If we have an enumeration, fill a dictionary.
    // The Python callpath handles this separately. The R callpath needs us
    // to do this. TODO: uniformize this at the callsites.
    if (column->has_enumeration() && use_enum) {
        auto enumvec = column->get_enumeration();

        ArrowSchema* dict_sch = new ArrowSchema;
        ArrowArray* dict_arr = new ArrowArray;

        dict_sch->format = (const char*)malloc(
            sizeof(char) * 2);  // mandatory, 'u' as 32bit indexing
        strcpy((char*)dict_sch->format, "u");
        dict_sch->name = nullptr;             // optional in dictionary
        dict_sch->metadata = nullptr;         // optional
        dict_sch->flags = 0;                  // optional
        dict_sch->n_children = 0;             // mandatory
        dict_sch->children = nullptr;         // optional
        dict_sch->dictionary = nullptr;       // optional
        dict_sch->release = &release_schema;  // mandatory
        dict_sch->private_data = nullptr;     // optional

        const int n_buf = 3;  // always variable here

        const int64_t n_vec = enumvec.size();
        dict_arr->length = n_vec;            // mandatory
        dict_arr->null_count = 0;            // mandatory
        dict_arr->offset = 0;                // mandatory
        dict_arr->n_buffers = n_buf;         // mandatory
        dict_arr->n_children = 0;            // mandatory
        dict_arr->buffers = nullptr;         // mandatory
        dict_arr->children = nullptr;        // optional
        dict_arr->dictionary = nullptr;      // optional
        dict_arr->release = &release_array;  // release from parent
        dict_arr->private_data = nullptr;    // optional here

        column->convert_enumeration();
        dict_arr->buffers = (const void**)malloc(sizeof(void*) * n_buf);
        dict_arr->buffers[0] = nullptr;  // validity: none here
        dict_arr->buffers[1] = column->enum_offsets().data();
        dict_arr->buffers[2] = column->enum_string().data();

        schema->dictionary = dict_sch;
        array->dictionary = dict_arr;
    }

    return std::pair(std::move(array), std::move(schema));
}

std::string_view ArrowAdapter::to_arrow_format(tiledb_datatype_t datatype) {
    switch (datatype) {
        case TILEDB_STRING_ASCII:
        case TILEDB_STRING_UTF8:
            return "U";  // large because TileDB uses 64bit offsets
        case TILEDB_CHAR:
        case TILEDB_BLOB:
            return "Z";  // large because TileDB uses 64bit offsets
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

}  // namespace tiledbsoma