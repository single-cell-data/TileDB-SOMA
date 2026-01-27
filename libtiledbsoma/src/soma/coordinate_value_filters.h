/**
 * @file   coordinate_value_filters.h
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

#ifndef TILEDBSOMA_COORDINATE_VALUE_FILTERS_H
#define TILEDBSOMA_COORDINATE_VALUE_FILTERS_H

#include <memory>
#include <optional>
#include <span>
#include <sstream>
#include <vector>

#include <tiledb/tiledb.h>
#include <tiledb/tiledb>

#include "../common/soma_column_selection.h"
#include "../soma/soma_column.h"
#include "../soma/soma_context.h"
#include "../tiledb_adapter/value_filter.h"
#include "../utils/common.h"

namespace tiledbsoma {
using namespace tiledb;

class CoordinateValueFilters {
   public:
    CoordinateValueFilters() = delete;
    CoordinateValueFilters(
        std::shared_ptr<Array> array,
        std::shared_ptr<SOMAContext> ctx,
        std::vector<std::shared_ptr<SOMAColumn>> index_columns,
        Domainish domain_kind);
    CoordinateValueFilters(CoordinateValueFilters&& other) = default;
    CoordinateValueFilters(const CoordinateValueFilters& other) = default;
    ~CoordinateValueFilters() = default;

    /**
     * @brief Add a constraint to a column.
     *
     * The query condition will select on elements of the TileDB array that include coordinates
     * within the selected range for the column (inclusive of end points). The slice must intersect
     * the domain.
     *
     * @tparam T The type of the column the query condition is applied to.
     * @param col_index The index of the column.
     * @param slice The coordinate selection to apply to the columns.
     */
    template <typename T>
    CoordinateValueFilters& add_slice(int64_t col_index, const SliceSelection<T>& selection) {
        auto col = index_columns_[col_index];
        if constexpr (std::is_same_v<T, std::string>) {
            validate_string_column(col);
        } else {
            try {
                tiledb::impl::type_check<T>(col->domain_type().value());
            } catch (const TileDBError& e) {
                throw std::invalid_argument("Invalid type.");  // TODO: Add error message before merging.
            }
            auto domain = col->domain_slot<T>(*ctx_, *array_, domain_kind_);
            if (!selection.has_overlap(domain)) {
                // Use sstream beacuse we don't want to include fmt directly in external header.
                std::stringstream ss;
                ss << "Non-overlapping slice [" << selection.start.value_or(domain.first) << ", "
                   << selection.stop.value_or(domain.second) << "] on column '" << col->name()
                   << "'. Slice must overlap the current column domain [" << domain.first << ", " << domain.second
                   << "].";
                throw std::out_of_range(ss.str());
            }
        }
        add_coordinate_query_condition(
            col_index, ValueFilter::create_from_slice<T>(*ctx_->tiledb_ctx(), col->name(), selection));

        return *this;
    }

    /**
     * @brief Add a constraint to a column.
     *
     * For a series of points: The query condition will select on elements of the TileDB array that
     * include coordinates that equal one of the provided elements for the specified attribute or
     * dimension. This method is not recommended for floating-point components. All points must be
     * inside the domain.
     *
     * @tparam T The type of the column the query condition is applied to.
     * @param col_index The index of the column.
     * @param points The coordinate selection to apply to the columns.
     */
    template <typename T>
    CoordinateValueFilters& add_points(int64_t col_index, PointSelection<T> selection) {
        auto col = index_columns_[col_index];
        if constexpr (std::is_same_v<T, std::string>) {
            validate_string_column(col);
        } else {
            try {
                tiledb::impl::type_check<T>(col->domain_type().value());
            } catch (const TileDBError& e) {
                throw std::invalid_argument("Invalid type.");  // TODO: Add error message before merging.
            }
            auto domain = col->domain_slot<T>(*ctx_, *array_, domain_kind_);
            if (!selection.is_subset(domain)) {
                // Use sstream beacuse we don't want to include fmt directly in external header.
                std::stringstream ss;
                ss << "Out-of-bounds coordinates found on column '" << col->name()
                   << "'. Coordinates must be inside column domain [" << domain.first << ", " << domain.second << "].";
                throw std::out_of_range(ss.str());
            }
        }
        return add_coordinate_query_condition(
            col_index, ValueFilter::create_from_points<T>(*ctx_->tiledb_ctx(), col->name(), selection));
    }

    /**
     * @brief Add a constraint to a column.
     *
     * For a slice constraint: the query condition will select on elements of the TileDB array that
     * include coordinates within the selected range for the column (inclusive of end points). The
     * slice must intersect the domain.
     *
     * For a series of points: The query condition will select on elements of the TileDB array that
     * include coordinates that equal one of the provided elements for the specified attribute or
     * dimension. This method is not recommended for floating-point components. All points must be
     * inside the domain.
     *
     * For monostate: no constraint is applied.
     *
     * @tparam T The type of the column the query condition is applied to.
     * @param col_index The index of the column.
     * @param selection The coordinate selection to apply to the columns.
     */
    template <typename T>
    CoordinateValueFilters& add_column_selection(int64_t col_index, SOMAColumnSelection<T> selection) {
        std::visit(
            [&](auto&& val) {
                using S = std::decay_t<decltype(val)>;
                if constexpr (std::is_same_v<S, SliceSelection<T>>) {
                    add_slice<T>(col_index, val);
                } else if constexpr (std::is_same_v<S, PointSelection<T>>) {
                    add_points<T>(col_index, val);
                }
                // Otherwise monostate: do nothing.
            },
            selection);
        return *this;
    }

    /**
     * Returns `true` if at least one query condition is set.
     */
    bool is_initialized() const;

    /**
     * Returns combined `ValueFilter` for the currently set coordinates.
     */
    ValueFilter combine() const;

   private:
    /**Add a query condition to a coordinate. */
    CoordinateValueFilters& add_coordinate_query_condition(int64_t index, ValueFilter&& qc);

    /**Throws an error if using a string on a non-string column. */
    void validate_string_column(std::shared_ptr<SOMAColumn> column) const;

    /** TileDB context required for creating query conditions. */
    std::shared_ptr<SOMAContext> ctx_;

    /** TileDB Array (only included because it's required to get domain from a column). */
    std::shared_ptr<Array> array_;

    /** Names of the dimensions of the underlying TileDB Array. */
    std::vector<std::shared_ptr<SOMAColumn>> index_columns_;

    /** The type of domain stored in the array. */
    Domainish domain_kind_;

    /** Query condition for coordinates. */
    std::vector<ValueFilter> coord_qc_;
};

}  // namespace tiledbsoma

#endif
