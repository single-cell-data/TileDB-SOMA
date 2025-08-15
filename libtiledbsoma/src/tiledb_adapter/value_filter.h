/**
 * @file   value_filter.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines helper functions and classes for using query conditions in SOMA.
 */

#ifndef TILEDBSOMA_VALUE_FILTER_H
#define TILEDBSOMA_VALUE_FILTER_H

#include <memory>
#include <optional>
#include <span>
#include <sstream>
#include <vector>

#include <tiledb/tiledb.h>
#include <tiledb/tiledb>

#include "../common/soma_column_selection.h"
#include "../utils/common.h"

namespace tiledbsoma {
using namespace tiledb;

class ValueFilter {
   public:
    ValueFilter() = default;
    ValueFilter(ValueFilter&& other) = default;
    ValueFilter(const ValueFilter& other) = default;
    ~ValueFilter() = default;

    ValueFilter& operator=(const ValueFilter& other) = default;
    ValueFilter& operator=(ValueFilter&& other) = default;

    /**
     * Create a ValueFilter from an initialized TileDB query condition.
     *
     * Important: The TileDB QueryCondtion must be initialized.
     *
     * @param qc The initialized TileDB query condition.
     */
    ValueFilter(const QueryCondition& qc);

    /**
     * Create a ValueFilter from an initialized TileDB query condition.
     *
     * Important: The TileDB QueryCondtion must be initialized.
     *
     * @param qc The initialized TileDB query condition.
     */
    ValueFilter(QueryCondition&& qc);

    /**
     * Create a query condition from a slice.
     *
     * @tparam The type of the element the query condition applies to.
     * @param column_name The name of the element the query condition applies to.
     * @param start_value The start of the range (inclusive).
     * @param stop_value The end of the range (inclusive).
     */
    template <typename T>
    static ValueFilter create_from_slice(
        const Context& ctx, const std::string& column_name, const SliceSelection<T>& slice) {
        if (slice.start.has_value() && slice.stop.has_value()) {
            return ValueFilter(
                QueryCondition::create(ctx, column_name, slice.start.value(), TILEDB_GE)
                    .combine(QueryCondition::create(ctx, column_name, slice.stop.value(), TILEDB_LE), TILEDB_AND));
        }
        if (slice.start.has_value()) {
            return ValueFilter(QueryCondition::create(ctx, column_name, slice.start.value(), TILEDB_GE));
        }
        if (slice.stop.has_value()) {
            return ValueFilter(QueryCondition::create(ctx, column_name, slice.stop.value(), TILEDB_LE));
        }
        return ValueFilter();
    }

    /**
     * Create a query condition from a series of points on a range.
     *
     * @tparam The type of the element the query condition applies to.
     * @param column_name The name of the element the query condition applies to.
     * @param values The values of the points
     */
    template <typename T>
    static ValueFilter create_from_points(
        const Context& ctx, const std::string& column_name, PointSelection<T> values) {
        if (values.points.empty()) {
            // Use sstream beacuse we don't want to include fmt directly in external header.
            std::stringstream ss;
            ss << "Cannot set coordinates on column '" << column_name << "'. No coordinates provided.";
            throw std::invalid_argument(ss.str());
        }
        if constexpr (std::is_same_v<T, std::string>) {
            // Using C API because C++ API only supports std::vector, not std::span.
            uint64_t data_size = 0;
            for (auto& val : values.points) {
                data_size += val.size();
            }
            std::vector<uint8_t> data(data_size);
            std::vector<uint64_t> offsets{};
            uint64_t curr_offset = 0;
            for (auto& val : values.points) {
                offsets.push_back(curr_offset);
                memcpy(data.data() + curr_offset, val.data(), val.size());
                curr_offset += val.size();
            }
            tiledb_query_condition_t* qc;
            ctx.handle_error(tiledb_query_condition_alloc_set_membership(
                ctx.ptr().get(),
                column_name.c_str(),
                data.data(),
                data.size(),
                offsets.data(),
                offsets.size() * sizeof(uint64_t),
                TILEDB_IN,
                &qc));
            return QueryCondition(ctx, qc);
        } else {
            // Using C API because C++ API only supports std::vector, not std::span.
            std::vector<uint64_t> offsets(values.points.size());
            for (size_t index = 0; index < values.points.size(); ++index) {
                offsets[index] = index * sizeof(T);
            }
            tiledb_query_condition_t* qc;
            ctx.handle_error(tiledb_query_condition_alloc_set_membership(
                ctx.ptr().get(),
                column_name.c_str(),
                values.points.data(),
                values.points.size() * sizeof(T),
                offsets.data(),
                offsets.size() * sizeof(uint64_t),
                TILEDB_IN,
                &qc));
            return QueryCondition(ctx, qc);
        }
    }

    /**
     * Return if the query condition is initialized.
     */
    inline bool is_initialized() const {
        return qc_.has_value();
    }

    /**
     * Return internal TileDB query condition.
     */
    inline const QueryCondition& query_condition() const {
        if (!qc_.has_value()) {
            throw TileDBSOMAError("Internal error: Cannot return query condition. Query condition is not initialized.");
        }
        return qc_.value();
    }

   private:
    /** Query condition for coordinates. */
    std::optional<QueryCondition> qc_;
};

}  // namespace tiledbsoma

#endif
