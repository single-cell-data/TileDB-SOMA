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
        std::shared_ptr<ColumnBuffer> column, bool use_enum = false) {
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

    /**
     * @brief Get Arrow format string from TileDB datatype.
     *
     * @param datatype TileDB datatype.
     * @return std::string_view Arrow format string.
     */
    static std::string_view to_arrow_format(tiledb_datatype_t datatype) {
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

    static ArraySchema arrow_schema_to_tiledb_schema(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        std::vector<std::string> index_column_names,
        std::vector<ArrowArray*> domain) {
        if (domain.size() != index_column_names.size()) {
            throw TileDBSOMAError(
                "if domain is specified, it must have the same length as "
                "index_column_names");
        }

        ArraySchema tdb_schema(*ctx, TILEDB_SPARSE);
        Domain tdb_dom(*ctx);

        for (size_t name_idx = 0; name_idx < index_column_names.size();
             ++name_idx) {
            for (int64_t schema_idx = 0; schema_idx < schema->n_children;
                 ++schema_idx) {
                auto child = schema->children[schema_idx];
                if (child->name == index_column_names[name_idx]) {
                    auto typeinfo = arrow_type_to_tiledb(child);
                    std::vector<std::byte> slot_domain;
                    std::vector<std::byte> tile_extent;

                    if (domain.size() != 0) {
                        auto data_buffer = domain[name_idx]->buffers[1];
                        auto dim_info = _schema_get_dim_from_buffer(
                            data_buffer, typeinfo.type);
                        slot_domain = std::get<0>(dim_info);
                        tile_extent = std::get<1>(dim_info);
                    } else {
                        auto dim_info = _schema_get_default_dim(child);
                        slot_domain = std::get<0>(dim_info);
                        tile_extent = std::get<1>(dim_info);
                    }

                    tdb_dom.add_dimension(Dimension::create(
                        *ctx,
                        child->name,
                        typeinfo.type,
                        slot_domain.data(),
                        tile_extent.data()));
                }
            }
        }
        tdb_schema.set_domain(tdb_dom);

        for (int64_t i = 0; i < schema->n_children; ++i) {
            auto child = schema->children[i];
            if (std::find(
                    index_column_names.begin(),
                    index_column_names.end(),
                    child->name) == index_column_names.end()) {
                auto typeinfo = arrow_type_to_tiledb(child);
                auto attr = Attribute(*ctx, child->name, typeinfo.type);
                tdb_schema.add_attribute(attr);
            }
        }

        // remember enums
        return tdb_schema;
    }

    struct TypeInfo {
        tiledb_datatype_t type;
        uint64_t elem_size;
        uint32_t cell_val_num;

        // is this represented as "Arrow large"
        bool arrow_large;
    };

    static TypeInfo arrow_type_to_tiledb(ArrowSchema* arw_schema) {
        auto fmt = std::string(arw_schema->format);
        bool large = false;
        if (fmt == "+l") {
            large = false;
            assert(arw_schema->n_children == 1);
            arw_schema = arw_schema->children[0];
        } else if (fmt == "+L") {
            large = true;
            assert(arw_schema->n_children == 1);
            arw_schema = arw_schema->children[0];
        }

        if (fmt == "i")
            return {TILEDB_INT32, 4, 1, large};
        else if (fmt == "l")
            return {TILEDB_INT64, 8, 1, large};
        else if (fmt == "f")
            return {TILEDB_FLOAT32, 4, 1, large};
        else if (fmt == "g")
            return {TILEDB_FLOAT64, 8, 1, large};
        else if (fmt == "b")
            return {TILEDB_BOOL, 1, 1, large};
        else if (fmt == "B")
            return {TILEDB_BLOB, 1, 1, large};
        else if (fmt == "c")
            return {TILEDB_INT8, 1, 1, large};
        else if (fmt == "C")
            return {TILEDB_UINT8, 1, 1, large};
        else if (fmt == "s")
            return {TILEDB_INT16, 2, 1, large};
        else if (fmt == "S")
            return {TILEDB_UINT16, 2, 1, large};
        else if (fmt == "I")
            return {TILEDB_UINT32, 4, 1, large};
        else if (fmt == "L")
            return {TILEDB_UINT64, 8, 1, large};
        // this is kind of a hack
        // technically 'tsn:' is timezone-specific, which we don't support
        // however, the blank (no suffix) base is interconvertible w/
        // np.datetime64
        else if (fmt == "tsn:")
            return {TILEDB_DATETIME_NS, 8, 1, large};
        else if (fmt == "z" || fmt == "Z")
            return {TILEDB_CHAR, 1, TILEDB_VAR_NUM, fmt == "Z"};
        else if (fmt == "u" || fmt == "U")
            return {TILEDB_STRING_UTF8, 1, TILEDB_VAR_NUM, fmt == "U"};
        else
            throw tiledb::TileDBError(
                "[TileDB-Arrow]: Unknown or unsupported Arrow format string '" +
                fmt + "'");
    };

    template <typename T>
    static std::vector<std::byte> __convert_domain_as_bytes(T min, T max) {
        std::byte* min_ptr = reinterpret_cast<std::byte*>(&min);
        std::byte* max_ptr = reinterpret_cast<std::byte*>(&max);
        std::vector<std::byte> slot_domain(min_ptr, min_ptr + sizeof(T));
        slot_domain.insert(slot_domain.end(), max_ptr, max_ptr + sizeof(T));
        return slot_domain;
    }

    template <typename T>
    static std::vector<std::byte> __convert_tile_as_bytes(T tile) {
        std::byte* tile_ptr = reinterpret_cast<std::byte*>(&tile);
        return std::vector<std::byte>(tile_ptr, tile_ptr + sizeof(T));
    }

    template <typename T>
    static std::tuple<std::vector<std::byte>, std::vector<std::byte>>
    __convert_dim_info_as_bytes(T min, T max) {
        T tile = std::min((T)2048, (T)(max - min + 1));
        return std::make_tuple(
            __convert_domain_as_bytes(min, max), __convert_tile_as_bytes(tile));
    }

    using dim_info_t =
        std::tuple<std::vector<std::byte>, std::vector<std::byte>>;

    template <typename T>
    static dim_info_t __convert_dim_info_as_bytes(T min, T max, T tile) {
        return std::make_tuple(
            __convert_domain_as_bytes(min, max), __convert_tile_as_bytes(tile));
    }

    static dim_info_t _schema_get_default_dim(ArrowSchema* arw_schema) {
        auto datatype = arrow_type_to_tiledb(arw_schema).type;
        switch (datatype) {
            case TILEDB_FLOAT32: {
                float min = std::numeric_limits<float>::min();
                float max = std::numeric_limits<float>::max();
                return __convert_dim_info_as_bytes(min, max, (float)2048);
            }
            case TILEDB_FLOAT64: {
                double min = std::numeric_limits<double>::min();
                double max = std::numeric_limits<double>::max();
                return __convert_dim_info_as_bytes(min, max, (double)2048);
            }
            case TILEDB_UINT8: {
                uint8_t min = std::numeric_limits<uint8_t>::min();
                uint8_t max = std::numeric_limits<uint8_t>::max();
                return __convert_dim_info_as_bytes(min, max, (uint8_t)64);
            }
            case TILEDB_INT8: {
                int8_t min = std::numeric_limits<int8_t>::min();
                int8_t max = std::numeric_limits<int8_t>::max();
                return __convert_dim_info_as_bytes(min, max, (int8_t)64);
            }
            case TILEDB_UINT16: {
                uint16_t min = std::numeric_limits<uint16_t>::min();
                uint16_t max = std::numeric_limits<uint16_t>::max();
                return __convert_dim_info_as_bytes(min, max, (uint16_t)2048);
            }
            case TILEDB_INT16: {
                int16_t min = std::numeric_limits<int16_t>::min();
                int16_t max = std::numeric_limits<int16_t>::max();
                return __convert_dim_info_as_bytes(min, max, (int16_t)2048);
            }
            case TILEDB_UINT32: {
                uint32_t min = std::numeric_limits<uint32_t>::min();
                uint32_t max = std::numeric_limits<uint32_t>::max();
                return __convert_dim_info_as_bytes(min, max, (uint32_t)2048);
            }
            case TILEDB_INT32: {
                int32_t min = std::numeric_limits<int32_t>::min();
                int32_t max = std::numeric_limits<int32_t>::max();
                return __convert_dim_info_as_bytes(min, max, (int32_t)2048);
            }
            case TILEDB_UINT64: {
                uint64_t min = std::numeric_limits<uint64_t>::min();
                uint64_t max = std::numeric_limits<uint64_t>::max();
                return __convert_dim_info_as_bytes(min, max, (uint64_t)2048);
            }
            case TILEDB_INT64:
            case TILEDB_DATETIME_YEAR:
            case TILEDB_DATETIME_MONTH:
            case TILEDB_DATETIME_WEEK:
            case TILEDB_DATETIME_DAY:
            case TILEDB_DATETIME_HR:
            case TILEDB_DATETIME_MIN:
            case TILEDB_DATETIME_SEC:
            case TILEDB_DATETIME_MS:
            case TILEDB_DATETIME_US:
            case TILEDB_DATETIME_NS:
            case TILEDB_DATETIME_PS:
            case TILEDB_DATETIME_FS:
            case TILEDB_DATETIME_AS: {
                int64_t min, max;
                if (strcmp(arw_schema->name, "soma_joinid") == 0) {
                    min = std::numeric_limits<uint32_t>::min();
                    max = std::numeric_limits<uint32_t>::max();
                } else {
                    min = std::numeric_limits<int64_t>::min();
                    max = std::numeric_limits<int64_t>::max();
                }
                return __convert_dim_info_as_bytes(min, max, (int64_t)2048);
            }
            default:
                throw tiledb::TileDBError(
                    "[TileDB-Arrow]: Unsupported TileDB type for dimension)");
        }
    };

    static dim_info_t _schema_get_dim_from_buffer(
        const void* buffer, tiledb_datatype_t datatype) {
        switch (datatype) {
            case TILEDB_FLOAT32: {
                float min = ((float*)buffer)[0];
                float max = ((float*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max);
            }
            case TILEDB_FLOAT64: {
                double min = ((double*)buffer)[0];
                double max = ((double*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max);
            }
            case TILEDB_UINT8: {
                uint8_t min = ((uint8_t*)buffer)[0];
                uint8_t max = ((uint8_t*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max, (uint8_t)64);
            }
            case TILEDB_INT8: {
                int8_t min = ((int8_t*)buffer)[0];
                int8_t max = ((int8_t*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max, (int8_t)64);
            }
            case TILEDB_UINT16: {
                uint16_t min = ((uint16_t*)buffer)[0];
                uint16_t max = ((uint16_t*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max);
            }
            case TILEDB_INT16: {
                int16_t min = ((int16_t*)buffer)[0];
                int16_t max = ((int16_t*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max);
            }
            case TILEDB_UINT32: {
                uint32_t min = ((uint32_t*)buffer)[0];
                uint32_t max = ((uint32_t*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max);
            }
            case TILEDB_INT32: {
                int32_t min = ((int32_t*)buffer)[0];
                int32_t max = ((int32_t*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max);
            }
            case TILEDB_UINT64: {
                uint64_t min = ((uint64_t*)buffer)[0];
                uint64_t max = ((uint64_t*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max);
            }
            case TILEDB_INT64:
            case TILEDB_DATETIME_YEAR:
            case TILEDB_DATETIME_MONTH:
            case TILEDB_DATETIME_WEEK:
            case TILEDB_DATETIME_DAY:
            case TILEDB_DATETIME_HR:
            case TILEDB_DATETIME_MIN:
            case TILEDB_DATETIME_SEC:
            case TILEDB_DATETIME_MS:
            case TILEDB_DATETIME_US:
            case TILEDB_DATETIME_NS:
            case TILEDB_DATETIME_PS:
            case TILEDB_DATETIME_FS:
            case TILEDB_DATETIME_AS: {
                int64_t min = ((int64_t*)buffer)[0];
                int64_t max = ((int64_t*)buffer)[1];
                return __convert_dim_info_as_bytes(min, max);
            }
            default:
                throw tiledb::TileDBError(
                    "[TileDB-Arrow]: Unsupported TileDB type for dimension)");
        }
    }
};
};  // namespace tiledbsoma

#endif
