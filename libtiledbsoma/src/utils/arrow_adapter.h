#ifndef ARROW_ADAPTER_H
#define ARROW_ADAPTER_H

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>

// https://arrow.apache.org/docs/format/CDataInterface.html
// https://arrow.apache.org/docs/format/Columnar.html#buffer-listing-for-each-layout
// https://arrow.apache.org/docs/format/CDataInterface.html#exporting-a-simple-int32-array

#ifndef ARROW_SCHEMA_AND_ARRAY_DEFINED
#include "carrow.h"
#endif
namespace tiledbsoma {

using namespace tiledb;

class ColumnBuffer;

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
    static void release_schema(struct ArrowSchema* schema);
    static void release_array(struct ArrowArray* array);

    /**
     * @brief Convert ColumnBuffer to an Arrow array.
     *
     * @return std::pair<std::unique_ptr<ArrowArray>,
     * std::unique_ptr<ArrowSchema>>
     */
    static std::pair<std::unique_ptr<ArrowArray>, std::unique_ptr<ArrowSchema>>
    to_arrow(std::shared_ptr<ColumnBuffer> column);

    static std::unique_ptr<ArrowSchema> arrow_schema_from_tiledb_array(
        std::shared_ptr<Context> ctx, std::shared_ptr<Array> tiledb_array);

    static std::unique_ptr<ArrowSchema> tiledb_schema_to_arrow_schema(
        std::shared_ptr<ArraySchema> tiledb_schema) {
        auto ndim = tiledb_schema->domain().ndim();
        auto nattr = tiledb_schema->attribute_num();

        std::unique_ptr<ArrowSchema>
            arrow_schema = std::make_unique<ArrowSchema>();
        arrow_schema->format = "+s";
        arrow_schema->n_children = ndim + nattr;
        arrow_schema->release = &release_schema;
        arrow_schema->children = (ArrowSchema**)malloc(
            sizeof(ArrowSchema*) * arrow_schema->n_children);

        ArrowSchema* child;

        for (uint32_t i = 0; i < ndim; ++i) {
            auto dim = tiledb_schema->domain().dimension(i);
            child = arrow_schema->children[i] = (ArrowSchema*)malloc(
                sizeof(ArrowSchema));
            child->format = to_arrow_format(dim.type()).data();
            child->name = strdup(dim.name().c_str());
            child->metadata = nullptr;
            child->flags = 0;
            child->n_children = 0;
            child->dictionary = nullptr;
            child->children = nullptr;
            child->release = &release_schema;
        }

        for (uint32_t i = 0; i < nattr; ++i) {
            auto attr = tiledb_schema->attribute(i);
            child = arrow_schema->children[ndim + i] = (ArrowSchema*)malloc(
                sizeof(ArrowSchema));
            child->format = to_arrow_format(attr.type()).data();
            child->name = strdup(attr.name().c_str());
            child->metadata = nullptr;
            child->flags = attr.nullable() ? ARROW_FLAG_NULLABLE : 0;
            child->n_children = 0;
            child->dictionary = nullptr;
            child->children = nullptr;
            child->release = &release_schema;
        }

        return arrow_schema;
    }

    /**
     * @brief Get Arrow format string from TileDB datatype.
     *
     * @param datatype TileDB datatype.
     * @return std::string_view Arrow format string.
     */
    static std::string_view to_arrow_format(
        tiledb_datatype_t datatype, bool use_large = true);

   private:
    static std::pair<const void*, std::size_t> _get_data_and_length(
        Enumeration& enmr, const void* dst);

    template <typename T>
    static const void* _fill_data_buffer(std::vector<T> src, const void* dst) {
        auto sz = src.size() * sizeof(T);
        dst = (const void*)malloc(sz);
        std::memcpy((void*)dst, src.data(), sz);
        return dst;
    }
};
};  // namespace tiledbsoma

#endif