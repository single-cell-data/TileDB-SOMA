#ifndef MANAGED_QUERY_H
#define MANAGED_QUERY_H

#include <tiledb/tiledb>

#include "tiledbsc/column_buffer.h"

namespace tiledbsc {

using namespace tiledb;

constexpr size_t TILEDBSC_DEFAULT_ALLOC = 524288;

class ManagedQuery {
   public:
    ManagedQuery(
        std::shared_ptr<Array> array,
        size_t initial_cells = TILEDBSC_DEFAULT_ALLOC);

    void select_columns(std::vector<std::string> names) {
        for (auto& name : names) {
            columns_.insert(name);
        }
    }

    template <typename T>
    void select_ranges(
        const std::string& dim, std::vector<std::pair<T, T>> ranges) {
        for (auto& [start, stop] : ranges) {
            query_->add_range(dim, start, stop);
        }
    }

    template <typename T>
    void select_points(const std::string& dim, std::vector<T> points) {
        for (auto& point : points) {
            query_->add_range(dim, point, point);
        }
    }

    size_t execute();

    template <typename T>
    std::span<T> data(const std::string& name) {
        return buffers_.at(name)->data<T>();
    }

    std::vector<std::string> strings(const std::string& name) {
        return buffers_.at(name)->strings();
    }

    std::string_view string_view(const std::string& name, uint64_t index) {
        return buffers_.at(name)->string_view(index);
    }

   private:
    std::shared_ptr<Array> array_;
    size_t initial_cells_;
    std::unique_ptr<Query> query_;
    std::set<std::string> columns_;
    std::unordered_map<std::string, std::shared_ptr<ColumnBuffer>> buffers_;
};

};  // namespace tiledbsc

#endif
