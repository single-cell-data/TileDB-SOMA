#include "tiledbsoma/soma_reader.h"
#include "tiledbsoma/logger_public.h"
#include "tiledbsoma/util.h"

namespace tiledbsoma {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMAReader> SOMAReader::open(
    std::string_view uri,
    std::string_view name,
    std::map<std::string, std::string> platform_config,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order) {
    return std::make_unique<SOMAReader>(
        uri,
        name,
        std::make_shared<Context>(Config(platform_config)),
        column_names,
        batch_size,
        result_order);
}

//===================================================================
//= public non-static
//===================================================================

SOMAReader::SOMAReader(
    std::string_view uri,
    std::string_view name,
    std::shared_ptr<Context> ctx,
    std::vector<std::string> column_names,
    std::string_view batch_size,
    std::string_view result_order)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(uri))
    , batch_size_(std::string(batch_size)) {
    // Validate parameters
    try {
        (void)result_order;
        LOG_DEBUG(fmt::format("[SOMAReader] opening array '{}'", uri_));
        auto array = std::make_shared<Array>(*ctx_, uri_, TILEDB_READ);
        mq_ = std::make_unique<ManagedQuery>(array, name);
        LOG_DEBUG(fmt::format("timestamp = {}", array->timestamp()));
    } catch (const std::exception& e) {
        throw TileDBSOMAError(
            fmt::format("Error opening array: {}\n  {}", uri_, e.what()));
    }

    if (!column_names.empty()) {
        mq_->select_columns(column_names);
    }
}

void SOMAReader::submit() {
    // Submit the query
    mq_->submit();
    first_read_next_ = true;
}

std::optional<std::shared_ptr<ArrayBuffers>> SOMAReader::read_next() {
    if (mq_->status() == Query::Status::UNINITIALIZED) {
        throw TileDBSOMAError(
            "[SOMAReader] submit must be called before read_next");
    }

    // Always return results from the first call to read_next()
    if (first_read_next_) {
        first_read_next_ = false;
        return mq_->results();
    }

    // If the query is complete, return `std::nullopt`.
    if (mq_->is_complete()) {
        return std::nullopt;
    }

    // Submit the query
    mq_->submit();

    // Return the results, possibly incomplete
    return mq_->results();
}

}  // namespace tiledbsoma
