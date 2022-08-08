#ifndef ARROW_ADAPTER_H
#define ARROW_ADAPTER_H

#include "carrow.h"
#include "tiledbsc/tiledbsc"

// https://arrow.apache.org/docs/format/CDataInterface.html
// https://arrow.apache.org/docs/format/Columnar.html#buffer-listing-for-each-layout
// https://arrow.apache.org/docs/format/CDataInterface.html#exporting-a-simple-int32-array

namespace tiledbsc {

class ArrowAdapter {
   public:
    static void release_schema(struct ArrowSchema* schema) {
        schema->release = nullptr;
    }

    static void release_array(struct ArrowArray* array) {
        free(array->buffers);
        array->release = nullptr;
    }

    /**
     * @brief Convert ColumnBuffer to an Arrow array.
     *
     * @return auto [Arrow array, Arrow schema]
     */
    static auto to_arrow(ColumnBuffer& column) {
        std::shared_ptr<ArrowSchema> schema = std::make_shared<ArrowSchema>();
        std::shared_ptr<ArrowArray> array = std::make_shared<ArrowArray>();

        schema->format = to_arrow_format(column.type()).data();  // mandatory
        schema->name = nullptr;                                  // optional
        schema->metadata = nullptr;                              // optional
        schema->flags = 0;                                       // optional
        schema->n_children = 0;                                  // mandatory
        schema->children = nullptr;                              // optional
        schema->dictionary = nullptr;                            // optional
        schema->release = &release_schema;                       // mandatory
        schema->private_data = nullptr;                          // optional

        int n_buffers = column.is_var() ? 3 : 2;

        array->length = column.size();    // mandatory
        array->null_count = 0;            // mandatory
        array->offset = 0;                // mandatory
        array->n_buffers = n_buffers;     // mandatory
        array->n_children = 0;            // mandatory
        array->buffers = nullptr;         // mandatory
        array->children = nullptr;        // optional
        array->dictionary = nullptr;      // optional
        array->release = &release_array;  // mandatory
        array->private_data = nullptr;    // mandatory

        array->buffers = (const void**)malloc(sizeof(void*) * n_buffers);
        assert(array->buffers != nullptr);
        array->buffers[0] = nullptr;  // validity
        if (n_buffers == 2) {
            array->buffers[1] = column.data<void*>().data();
        } else if (n_buffers == 3) {
            array->buffers[1] = column.offsets().data();
            array->buffers[2] = column.data<void*>().data();
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
                return "C";  // TILEDB_BOOL is 8bit but arrow BOOL is 1bit
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
        throw TileDBSCError(fmt::format(
            "ArrowAdapter: Unsupported TileDB datatype: {} ",
            tiledb::impl::type_to_str(datatype)));
    }
};

};  // namespace tiledbsc

#endif
