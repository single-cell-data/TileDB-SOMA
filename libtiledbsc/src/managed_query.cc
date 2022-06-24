#include "tiledbsc/managed_query.h"
#include "tiledbsc/common.h"
#include "tiledbsc/logger_public.h"

namespace tiledbsc {

ManagedQuery::ManagedQuery(std::shared_ptr<Array> array, size_t initial_cells)
    : array_(array)
    , initial_cells_(initial_cells) {
    query_ = std::make_unique<Query>(array->schema().context(), *array);

    if (array->schema().array_type() == TILEDB_SPARSE) {
        query_->set_layout(TILEDB_UNORDERED);
    } else {
        query_->set_layout(TILEDB_ROW_MAJOR);
    }
}

size_t ManagedQuery::execute() {
    auto status = query_->query_status();

    // Query is complete, return 0 cells read
    if (status == Query::Status::COMPLETE) {
        return 0;
    }

    // Query is uninitialized, allocate and attach buffers
    if (status == Query::Status::UNINITIALIZED) {
        // If no columns were selected, select all columns
        if (!columns_.size()) {
            for (const auto& dim : array_->schema().domain().dimensions()) {
                columns_.insert(dim.name());
            }
            for (const auto& [name, attr] : array_->schema().attributes()) {
                (void)attr;
                columns_.insert(name);
            }
        }

        // Allocate and attach buffers
        for (auto& name : columns_) {
            LOG_DEBUG(fmt::format("Adding buffer for column '{}'", name));
            buffers_.emplace(
                name, ColumnBuffer::create(array_, name, initial_cells_));
            buffers_[name]->attach(*query_);
        }
    }

    // Submit query
    query_->submit();
    status = query_->query_status();
    LOG_DEBUG(fmt::format("Query status = {}", (int)status));

    // Update ColumnBuffer size to match query results
    size_t num_cells = 0;
    for (auto& [name, buffer] : buffers_) {
        num_cells = buffer->update_size(*query_);
        LOG_DEBUG(fmt::format("Buffer {} cells={}", name, num_cells));
    }

    // TODO: retry the query with larger buffers
    if (status == Query::Status::INCOMPLETE && !num_cells) {
        throw TileDBSCError(fmt::format(
            "[ManagedQuery] Buffers are too small: {} cells", initial_cells_));
    }

    return num_cells;
}

};  // namespace tiledbsc
