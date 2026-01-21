#include "exporter.h"

#include "arrow_buffer.h"

#include "../logging/impl/logger.h"
#include "../logging/logger.h"

namespace tiledbsoma::common::arrow {
ArrowTable column_buffer_to_arrow(
    ColumnBuffer* column_buffer,
    const std::unordered_map<std::string, std::shared_future<std::shared_ptr<ArrowBuffer>>>& enumerations,
    bool downcast_dict_of_large_var) {
    auto [array, schema] = make_empty_arrow_table(column_buffer->name(), to_arrow_format(column_buffer->type()), 0);

    // this will be 3 for char vecs and 2 for enumerations
    array->n_buffers = column_buffer->is_var() ? 3 : 2;

    // Create an ArrowBuffer to manage the lifetime of `column`.
    // - `arrow_buffer` holds shared_ptr to `column`, increments
    //   the use count and keeps the ColumnBuffer data alive.
    // - When the arrow array is released, `array->release()` is
    //   called with `arrow_buffer` in `private_data`.
    //   `arrow_buffer` is deleted, which decrements the the
    //   `column` use count. When the `column` use count reaches
    //   0, the ColumnBuffer data will be deleted.
    std::shared_ptr<ArrowBuffer> buffer = std::make_shared<ArrowBuffer>(
        column_buffer->export_buffers(), column_buffer->name());

    array->length = buffer->storage()->length();

    if (column_buffer->is_nullable()) {
        schema->flags |= ARROW_FLAG_NULLABLE;
        array->null_count = buffer->storage()->null_count();
    } else {
        schema->flags &= ~ARROW_FLAG_NULLABLE;
    }

    logging::LOG_TRACE(
        fmt::format(
            "[column_buffer_to_arrow] column type {} name {} nbuf {} nullable {}",
            to_arrow_format(column_buffer->type()).data(),
            column_buffer->name().data(),
            array->n_buffers,
            column_buffer->is_nullable()));

    reinterpret_cast<PrivateArrowBuffer*>(array->private_data)->buffer_ = buffer;

    logging::LOG_TRACE(fmt::format("[column_buffer_to_arrow] create array name='{}'", column_buffer->name()));

    assert(array->buffers != nullptr);
    array->buffers[0] = nullptr;  // validity addressed below
    array->buffers[array->n_buffers - 1] = buffer->storage()->data().data();
    if (column_buffer->is_var()) {
        array->buffers[1] = buffer->storage()->offsets().data();
    }

    if (column_buffer->is_nullable()) {
        array->buffers[0] = buffer->storage()->validity().data();
    }

    auto enumeration = column_buffer->get_enumeration();
    if (enumeration.has_value()) {
        if (column_buffer->is_ordered()) {
            schema->flags |= ARROW_FLAG_DICTIONARY_ORDERED;
        }

        std::shared_ptr<ArrowBuffer> enumeration_buffer = enumerations.at(enumeration->name()).get();

        auto [dictionary_array, dictionary_schema] = make_empty_arrow_table(
            enumeration->name(), to_arrow_format(enumeration->type(), !downcast_dict_of_large_var), 0);

        bool is_var_enum = enumeration->type() == TILEDB_STRING_ASCII || enumeration->type() == TILEDB_STRING_UTF8 ||
                           enumeration->type() == TILEDB_CHAR || enumeration->type() == TILEDB_BLOB;
        dictionary_array->n_buffers = is_var_enum ? 3 : 2;
        assert(dictionary_array->buffers != nullptr);

        dictionary_array->buffers[0] = nullptr;
        dictionary_array->buffers[dictionary_array->n_buffers - 1] = enumeration_buffer->storage()->data().data();
        if (is_var_enum) {
            dictionary_array->buffers[1] = enumeration_buffer->storage()->offsets().data();
        }

        dictionary_array->length = enumeration_buffer->storage()->length();
        reinterpret_cast<PrivateArrowBuffer*>(dictionary_array->private_data)->buffer_ = enumeration_buffer;

        schema->dictionary = dictionary_schema.release();
        array->dictionary = dictionary_array.release();
    }

    return std::pair(std::move(array), std::move(schema));
}
}  // namespace tiledbsoma::common::arrow