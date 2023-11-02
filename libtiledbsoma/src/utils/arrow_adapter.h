#ifndef ARROW_ADAPTER_H
#define ARROW_ADAPTER_H

#include <tiledb/tiledb>
#include "../soma/column_buffer.h"
#include "../utils/logger.h"
#ifndef ARROW_SCHEMA_AND_ARRAY_DEFINED
#include "carrow.h"
#endif

// https://arrow.apache.org/docs/format/CDataInterface.html
// https://arrow.apache.org/docs/format/Columnar.html#buffer-listing-for-each-layout
// https://arrow.apache.org/docs/format/CDataInterface.html#exporting-a-simple-int32-array

namespace tiledbsoma {

/**
 * @brief The ArrowBuffer holds a shared pointer to a ColumnBuffer, which
 * manages the lifetime of a ColumnBuffer used to back an Arrow array.
 *
 * The ArrowArray.release callback will delete the ArrowBuffer, and
 * automatically decrement the use count of the ColumnBuffer's shared pointer.
 *
 */
struct ArrowBuffer {
    ArrowBuffer(std::shared_ptr<ColumnBuffer> buffer)
        : buffer_(buffer){};

    std::shared_ptr<ColumnBuffer> buffer_;
};

class ArrowAdapter {
   public:
    static void release_schema(struct ArrowSchema* schema) {
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

    static void release_array(struct ArrowArray* array) {
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

    /**
     * @brief Convert ColumnBuffer to an Arrow array.
     *
     * @return auto [Arrow array, Arrow schema]
     */
    static auto to_arrow(
        std::shared_ptr<ColumnBuffer> column, bool use_enum = true) {
        std::unique_ptr<ArrowSchema> schema = std::make_unique<ArrowSchema>();
        std::unique_ptr<ArrowArray> array = std::make_unique<ArrowArray>();

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
        // - `arrow_buffer` holds a shared_ptr to `column`, which increments
        //   the use count and keeps the ColumnBuffer data alive.
        // - When the arrow array is released, `array->release()` is called with
        //   `arrow_buffer` in `private_data`. `arrow_buffer` is deleted, which
        //   decrements the the `column` use count. When the `column` use count
        //   reaches 0, the ColumnBuffer data will be deleted.
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

        array->buffers = (const void**)malloc(sizeof(void*) * n_buffers);
        assert(array->buffers != nullptr);
        array->buffers[0] = nullptr;  // validity
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

        if (column->has_enumeration() && use_enum) {
            ArrowSchema* dict_sch = new ArrowSchema;
            ArrowArray* dict_arr = new ArrowArray;

            auto enmr = column->get_enumeration_info();
            dict_sch->format = strdup(
                to_arrow_format(enmr->type(), false).data());
            dict_sch->name = strdup(enmr->name().c_str());
            dict_sch->metadata = nullptr;
            dict_sch->flags = 0;
            dict_sch->n_children = 0;
            dict_sch->children = nullptr;
            dict_sch->dictionary = nullptr;
            dict_sch->release = &release_schema;
            dict_sch->private_data = nullptr;

            const int n_buf = strcmp(dict_sch->format, "u") == 0 ? 3 : 2;
            dict_arr->null_count = 0;
            dict_arr->offset = 0;
            dict_arr->n_buffers = n_buf;
            dict_arr->n_children = 0;
            dict_arr->buffers = nullptr;
            dict_arr->children = nullptr;
            dict_arr->dictionary = nullptr;
            dict_arr->release = &release_array;
            dict_arr->private_data = nullptr;

            dict_arr->buffers = (const void**)malloc(sizeof(void*) * n_buf);
            dict_arr->buffers[0] = nullptr;  // validity: none here

            // TODO string types currently get the data and offset buffers from
            // ColumnBuffer::enum_offsets and ColumnBuffer::enum_string which is
            // retrieved via ColumnBuffer::convert_enumeration. This may be
            // refactored to all use ColumnBuffer::get_enumeration_info. Note
            // that ColumnBuffer::has_enumeration may also be removed in a
            // future refactor as ColumnBuffer::get_enumeration_info returns
            // std::optional where std::nullopt indicates the column does not
            // contain enumerated values.
            if (enmr->type() == TILEDB_STRING_ASCII or
                enmr->type() == TILEDB_STRING_UTF8) {
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

    static std::pair<const void*, std::size_t> _get_data_and_length(
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
                dst = (const void*)malloc(sz);
                std::memcpy((void*)dst, &src, sz);

                return std::pair(dst, data.size());
            }
            case TILEDB_INT8: {
                auto data = enmr.as_vector<int8_t>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_UINT8: {
                auto data = enmr.as_vector<uint8_t>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_INT16: {
                auto data = enmr.as_vector<int16_t>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_UINT16: {
                auto data = enmr.as_vector<uint16_t>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_INT32: {
                auto data = enmr.as_vector<int32_t>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_UINT32: {
                auto data = enmr.as_vector<uint32_t>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_INT64: {
                auto data = enmr.as_vector<int64_t>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_UINT64: {
                auto data = enmr.as_vector<uint64_t>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_FLOAT32: {
                auto data = enmr.as_vector<float>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            case TILEDB_FLOAT64: {
                auto data = enmr.as_vector<double>();
                return std::pair(_fill_data_buffer(data, dst), data.size());
            }
            default:
                throw TileDBSOMAError(fmt::format(
                    "ArrowAdapter: Unsupported TileDB dict datatype: {} ",
                    tiledb::impl::type_to_str(enmr.type())));
        }
    }

    template <typename T>
    static const void* _fill_data_buffer(std::vector<T> src, const void* dst) {
        auto sz = src.size() * sizeof(T);
        dst = (const void*)malloc(sz);
        std::memcpy((void*)dst, src.data(), sz);
        return dst;
    }

    /**
     * @brief Get Arrow format string from TileDB datatype.
     *
     * @param datatype TileDB datatype.
     * @return std::string_view Arrow format string.
     */
    static std::string_view to_arrow_format(
        tiledb_datatype_t datatype, bool use_large = true) {
        switch (datatype) {
            case TILEDB_STRING_ASCII:
            case TILEDB_STRING_UTF8:
                return use_large ?
                           "U" :
                           "u";  // large because TileDB uses 64bit offsets
            case TILEDB_CHAR:
            case TILEDB_BLOB:
                return use_large ?
                           "Z" :
                           "z";  // large because TileDB uses 64bit offsets
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
};  // namespace tiledbsoma
};  // namespace tiledbsoma

#endif