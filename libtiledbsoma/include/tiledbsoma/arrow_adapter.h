#ifndef ARROW_ADAPTER_H
#define ARROW_ADAPTER_H

#include <tiledbsoma/tiledbsoma>
#include "carrow.h"

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
        array->release = nullptr;
    }

    /**
     * @brief Convert ColumnBuffer to an Arrow array.
     *
     * @return auto [Arrow array, Arrow schema]
     */
    static auto to_arrow(std::shared_ptr<ColumnBuffer> column) {
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
};

};  // namespace tiledbsoma

#endif
