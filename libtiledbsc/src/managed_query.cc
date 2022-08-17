#include "tiledbsc/managed_query.h"
#include "tiledbsc/common.h"
#include "tiledbsc/logger_public.h"

namespace tiledbsc {

using namespace tiledb;

//===================================================================
//= public non-static
//===================================================================

ManagedQuery::ManagedQuery(std::shared_ptr<Array> array, std::string name)
    : array_(array)
    , name_(name)
    , schema_(array->schema()) {
    query_ = std::make_unique<Query>(schema_.context(), *array);
    subarray_ = std::make_unique<Subarray>(schema_.context(), *array);

    if (array->schema().array_type() == TILEDB_SPARSE) {
        query_->set_layout(TILEDB_UNORDERED);
    } else {
        query_->set_layout(TILEDB_ROW_MAJOR);
    }
}

void ManagedQuery::select_columns(
    std::vector<std::string> names, bool if_not_empty) {
    // Return if we are selecting all columns (columns_ is empty) and we want to
    // continue selecting all columns (if_not_empty == true).
    if (if_not_empty && columns_.empty()) {
        return;
    }

    for (auto& name : names) {
        // Name is not an attribute or dimension.
        if (!schema_.has_attribute(name) &&
            !schema_.domain().has_dimension(name)) {
            LOG_DEBUG(fmt::format(
                "[ManagedQuery] [{}] Invalid column selected: {}",
                name_,
                name));
            invalid_columns_selected_ = true;
        } else {
            columns_.insert(name);
        }
    }
}

size_t ManagedQuery::submit() {
    auto status = query_->query_status();

    // Query is complete, return 0 cells read
    if (status == Query::Status::COMPLETE) {
        return 0;
    }

    // Query is uninitialized, allocate and attach buffers
    if (status == Query::Status::UNINITIALIZED) {
        // Set the subarray for range slicing
        query_->set_subarray(*subarray_);

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
            LOG_DEBUG(fmt::format(
                "[ManagedQuery] [{}] Adding buffer for column '{}'",
                name_,
                name));
            buffers_.emplace(name, ColumnBuffer::create(array_, name));
            buffers_[name]->attach(*query_);
        }
    }

    // Submit query
    LOG_DEBUG(fmt::format("[ManagedQuery] [{}] Submit query", name_));
    query_->submit();
    status = query_->query_status();
    LOG_DEBUG(fmt::format(
        "[ManagedQuery] [{}] Query status = {}", name_, (int)status));

    // If the query was ever incomplete, the result buffers contents are not
    // complete.
    if (status == Query::Status::INCOMPLETE) {
        results_complete_ = false;
    }

    // Update ColumnBuffer size to match query results
    size_t num_cells = 0;
    for (auto& [name, buffer] : buffers_) {
        num_cells = buffer->update_size(*query_);
        LOG_DEBUG(fmt::format(
            "[ManagedQuery] [{}] Buffer {} cells={}", name_, name, num_cells));
    }
    total_num_cells_ += num_cells;

    // TODO: retry the query with larger buffers
    if (status == Query::Status::INCOMPLETE && !num_cells) {
        throw TileDBSCError(
            fmt::format("[ManagedQuery] [{}] Buffers are too small.", name_));
    }

    return num_cells;
}

};  // namespace tiledbsc
