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

struct TypeInfo {
    tiledb_datatype_t type;
    uint64_t elem_size;
    uint32_t cell_val_num;

    // is this represented as "Arrow large"
    bool arrow_large;
};

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

    static ArraySchema arrow_schema_to_tiledb_schema(
        std::shared_ptr<Context> ctx,
        ArrowSchema& schema,
        std::vector<std::string> index_column_names,
        ArrowArray& domains,
        ArrowArray& extents) {
        if (domains.n_children != (int64_t)index_column_names.size()) {
            throw TileDBSOMAError(
                "if domain is specified, it must have the same length as "
                "index_column_names");
        }

        ArraySchema tdb_schema(*ctx, TILEDB_SPARSE);
        Domain tdb_dom(*ctx);

        auto platform_config = ctx->config();

        for (size_t col_idx = 0; col_idx < index_column_names.size();
             ++col_idx) {
            for (int64_t schema_idx = 0; schema_idx < schema.n_children;
                 ++schema_idx) {
                auto child = schema.children[schema_idx];
                auto typeinfo = arrow_type_to_tiledb(child);
                if (child->name == index_column_names[col_idx]) {
                    auto dim = Dimension::create(
                        *ctx,
                        child->name,
                        typeinfo.type,
                        domains.children[col_idx]->buffers[1],
                        extents.children[col_idx]->buffers[1]);

                    Filter filter(*ctx, TILEDB_FILTER_ZSTD);
                    if (platform_config.contains("dataframe_dim_zstd_level")) {
                        auto level = std::stoi(
                            platform_config.get("dataframe_dim_zstd_level"));
                        filter.set_option(TILEDB_COMPRESSION_LEVEL, level);
                    } else {
                        filter.set_option(TILEDB_COMPRESSION_LEVEL, 3);
                    }

                    FilterList filter_list(*ctx);
                    filter_list.add_filter(filter);
                    dim.set_filter_list(filter_list);

                    tdb_dom.add_dimension(dim);
                } else {
                    auto attr = Attribute(*ctx, child->name, typeinfo.type);
                    // std::cout << "child->flags: " << child->flags <<
                    // std::endl; std::cout << "ARROW_FLAG_NULLABLE: " <<
                    // ARROW_FLAG_NULLABLE
                    //           << std::endl;
                    // std::cout << child->flags | ARROW_FLAG_NULLABLE
                    //                                 << std::endl;
                    // if (child->flags | ARROW_FLAG_NULLABLE) {
                    //     attr.set_nullable(true);
                    // }
                    tdb_schema.add_attribute(attr);
                }
            }
        }

        tdb_schema.set_domain(tdb_dom);

        if (platform_config.contains("offsets_filters")) {
            auto filter_list = _get_filterlist_from_config_option(
                ctx, platform_config.get("offsets_filters"));
            tdb_schema.set_offsets_filter_list(filter_list);
        } else {
            FilterList filter_list(*ctx);
            filter_list.add_filter(Filter(*ctx, TILEDB_FILTER_DOUBLE_DELTA));
            filter_list.add_filter(
                Filter(*ctx, TILEDB_FILTER_BIT_WIDTH_REDUCTION));
            filter_list.add_filter(Filter(*ctx, TILEDB_FILTER_ZSTD));
            tdb_schema.set_offsets_filter_list(filter_list);
        }

        if (platform_config.contains("validity_filters")) {
            auto filter_list = _get_filterlist_from_config_option(
                ctx, platform_config.get("validity_filters"));
            tdb_schema.set_validity_filter_list(filter_list);
        }

        if (platform_config.contains("tile_order")) {
            auto order = _get_order_from_config_option(
                platform_config.get("tile_order"));
            tdb_schema.set_tile_order(order);
        }

        if (platform_config.contains("cell_order")) {
            auto order = _get_order_from_config_option(
                platform_config.get("cell_order"));
            tdb_schema.set_cell_order(order);
        }

        if (platform_config.contains("capacity")) {
            tdb_schema.set_capacity(std::stoi(platform_config.get("capacity")));
        }

        if (platform_config.contains("allows_duplicates")) {
            bool allows_duplicates;
            std::istringstream(platform_config.get("allows_duplicates")) >>
                std::boolalpha >> allows_duplicates;
            tdb_schema.set_allows_dups(allows_duplicates);
        }

        tdb_schema.check();

        return tdb_schema;
    }

    static FilterList _get_filterlist_from_config_option(
        std::shared_ptr<Context> ctx, std::string filters) {
        std::map<std::string, tiledb_filter_type_t> token_to_filter = {
            {"GzipFilter", TILEDB_FILTER_GZIP},
            {"ZstdFilter", TILEDB_FILTER_ZSTD},
            {"LZ4Filter", TILEDB_FILTER_LZ4},
            {"Bzip2Filter", TILEDB_FILTER_BZIP2},
            {"RleFilter", TILEDB_FILTER_RLE},
            {"DeltaFilter", TILEDB_FILTER_DELTA},
            {"DoubleDeltaFilter", TILEDB_FILTER_DOUBLE_DELTA},
            {"BitWidthReductionFilter", TILEDB_FILTER_BIT_WIDTH_REDUCTION},
            {"BitShuffleFilter", TILEDB_FILTER_BITSHUFFLE},
            {"ByteShuffleFilter", TILEDB_FILTER_BYTESHUFFLE},
            {"PositiveDeltaFilter", TILEDB_FILTER_POSITIVE_DELTA},
            {"ChecksumMD5Filter", TILEDB_FILTER_CHECKSUM_MD5},
            {"ChecksumSHA256Filter", TILEDB_FILTER_CHECKSUM_SHA256},
            {"DictionaryFilter", TILEDB_FILTER_DICTIONARY},
            {"FloatScaleFilter", TILEDB_FILTER_SCALE_FLOAT},
            {"XORFilter", TILEDB_FILTER_XOR},
            {"WebpFilter", TILEDB_FILTER_WEBP},
            {"NoOpFilter", TILEDB_FILTER_NONE}};

        FilterList filter_list(*ctx);
        std::istringstream ss(filters);
        std::string filter_tokenized;

        while (std::getline(ss, filter_tokenized, ',')) {
            filter_list.add_filter(
                Filter(*ctx, token_to_filter[filter_tokenized]));
        }
        return filter_list;
    }

    static tiledb_layout_t _get_order_from_config_option(std::string order) {
        std::map<std::string, tiledb_layout_t> token_to_order = {
            {"row-major", TILEDB_ROW_MAJOR},
            {"R", TILEDB_ROW_MAJOR},
            {"col-major", TILEDB_COL_MAJOR},
            {"C", TILEDB_COL_MAJOR},
            {"hilbert", TILEDB_HILBERT}};
        return token_to_order[order];
    }
};
};  // namespace tiledbsoma

#endif
