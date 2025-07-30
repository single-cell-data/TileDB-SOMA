/**
 * @file   soma_query_condition.h
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

#ifndef TILEDBSOMA_TILEDB_ADAPTER_H
#define TILEDBSOMA_TILEDB_ADAPTER_H

#include <memory>
#include <optional>
#include <span>
#include <sstream>
#include <vector>

#include <tiledb/tiledb.h>
#include <tiledb/tiledb>

#include "../common/soma_column_selection.h"
#include "../soma/soma_context.h"
#include "../utils/common.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAValueFilter {
   public:
    SOMAValueFilter() = default;
    SOMAValueFilter(SOMAValueFilter&& other) = default;
    SOMAValueFilter(const SOMAValueFilter& other) = default;
    ~SOMAValueFilter() = default;

    SOMAValueFilter& operator=(const SOMAValueFilter& other) = default;
    SOMAValueFilter& operator=(SOMAValueFilter&& other) = default;

    /**
     * Create a SOMAValueFilter from an initialized TileDB query condition.
     *
     * Important: The TileDB QueryCondtion must be initialized.
     *
     * @param qc The initialized TileDB query condition.
     */
    SOMAValueFilter(const QueryCondition& qc);

    /**
     * Create a SOMAValueFilter from an initialized TileDB query condition.
     *
     * Important: The TileDB QueryCondtion must be initialized.
     *
     * @param qc The initialized TileDB query condition.
     */
    SOMAValueFilter(QueryCondition&& qc);

    /**
     * Create a query condition from a slice.
     *
     * @tparam The type of the element the query condition applies to.
     * @param column_name The name of the element the query condition applies to.
     * @param start_value The start of the range (inclusive).
     * @param stop_value The end of the range (inclusive).
     */
    template <typename T, typename = std::enable_if<!std::is_same_v<T, std::string>>>
    static SOMAValueFilter create_from_range(
        const Context& ctx, const std::string& column_name, T start_value, T stop_value) {
        if (stop_value < start_value) {
            // Use sstream beacuse we don't want to include fmt directly in external header.
            std::stringstream ss;
            ss << "Cannot set range [" << start_value << ", " << stop_value << "] on column '" << column_name
               << "'. Invalid range: the final value must be greater "
                  "than or equal to the starting value.";
            throw std::invalid_argument(ss.str());
        }
        return SOMAValueFilter(
            QueryCondition::create<T>(ctx, column_name, start_value, TILEDB_GE)
                .combine(QueryCondition::create<T>(ctx, column_name, stop_value, TILEDB_LE), TILEDB_AND));
    }

    /**
     * Create a query condition from a slice.
     *
     * @param column_name The name of the element the query condition applies to.
     * @param start_value The start of the range (inclusive).
     * @param stop_value The end of the range (inclusive).
     */
    static SOMAValueFilter create_from_range(
        const Context& ctx,
        const std::string& column_name,
        const std::string& start_value,
        const std::string& stop_value);

    /**
     * Create a query condition from a series of points on a range.
     *
     * @tparam The type of the element the query condition applies to.
     * @param column_name The name of the element the query condition applies to.
     * @param values The values of the points
     */
    template <typename T, typename = std::enable_if<!std::is_same_v<T, std::string>>>
    static SOMAValueFilter create_from_points(const Context& ctx, const std::string& column_name, std::span<T> values) {
        if (values.empty()) {
            // Use sstream beacuse we don't want to include fmt directly in external header.
            std::stringstream ss;
            ss << "Cannot set coordinates on column '" << column_name << "'. No coordinates provided.";
            throw std::invalid_argument(ss.str());
        }
        // Using C API because C++ API only supports std::vector, not std::span.
        std::vector<uint64_t> offsets(values.size());
        for (size_t index = 0; index < values.size(); ++index) {
            offsets[index] = index * sizeof(T);
        }
        tiledb_query_condition_t* qc;
        ctx.handle_error(tiledb_query_condition_alloc_set_membership(
            ctx.ptr().get(),
            column_name.c_str(),
            values.data(),
            values.size() * sizeof(T),
            offsets.data(),
            offsets.size() * sizeof(uint64_t),
            TILEDB_IN,
            &qc));
        return QueryCondition(ctx, qc);
    }

    static SOMAValueFilter create_from_points(
        const Context& ctx, const std::string& column_name, std::span<std::string> values);
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

class SOMAIndexValueFilter {
   public:
    SOMAIndexValueFilter() = delete;
    SOMAIndexValueFilter(const SOMAContext& ctx, const std::vector<std::string>& dim_names);

    SOMAIndexValueFilter(SOMAIndexValueFilter&& other) = default;
    SOMAIndexValueFilter(const SOMAIndexValueFilter& other) = default;
    ~SOMAIndexValueFilter() = default;

    /**
     * @brief Add a constraint to a column.
     *
     * For a slice constraint: the query condition will select on elements of the TileDB array that 
     * include coordinates within the selected range for the column (inclusive of end points). If a
     * domain is provided, the slice must intersect the domain.
     *
     * For a series of points: The query condition will select on elements of the TileDB array that 
     * include coordinates that equal one of the provided elements for the specified attribute or 
     * dimension. This method is not recommended for floating-point components. If a domain is provided,
     * all points must be inside the domain.
     *
     * For monostate: no constraint is applied.
     *
     * @tparam T The type of the column the query condition is applied to.
     * @param col_index The index of the column 
     * @param selection The coordinate selection to apply to the columns.
     * @param domain Optional domain of the column the selection is applied to.
     */
    template <typename T, typename = std::enable_if<!std::is_same_v<T, std::string>>>
    SOMAIndexValueFilter& add_column_selection(
        int64_t col_index, SOMAColumnSelection<T> selection, const std::pair<T, T>& domain) {
        std::visit(
            [&](auto&& val) {
                using S = std::decay_t<decltype(val)>;
                if constexpr (std::is_same_v<S, SOMASliceSelection<T>>) {
                    if (!val.has_overlap(domain)) {
                        // Use sstream beacuse we don't want to include fmt directly in external header.
                        std::stringstream ss;
                        ss << "Non-overlapping slice [" << val.start << ", " << val.stop << "] on column '"
                           << dim_names_[col_index] << "'. Slice must overlap the current column domain ["
                           << domain.first << ", " << domain.second << "].";
                        throw std::out_of_range(ss.str());
                    }
                    add_coordinate_query_condition(
                        col_index,
                        SOMAValueFilter::create_from_range<T>(*ctx_, dim_names_[col_index], val.start, val.stop));

                } else if constexpr (std::is_same_v<S, SOMAPointSelection<T>>) {
                    if (!val.is_subset(domain)) {
                        // Use sstream beacuse we don't want to include fmt directly in external header.
                        std::stringstream ss;
                        ss << "Out-of-bounds coordinates found on column '" << dim_names_[col_index]
                           << "'. Coordinates must be inside column domain [" << domain.first << ", " << domain.second
                           << "].";
                        throw std::out_of_range(ss.str());
                    }
                    add_coordinate_query_condition(
                        col_index, SOMAValueFilter::create_from_points<T>(*ctx_, dim_names_[col_index], val.points));
                }
                // Otherwise monostate: do nothing.
            },
            selection);
        return *this;
    }

    /**
     * @brief Add a constraint to a string column.
     *
     * For a slice constraint: the query condition will select on elements of the TileDB array that 
     * include coordinates within the selected range for the column (inclusive of end points).
     *
     * For a series of points: The query condition will select on elements of the TileDB array that 
     * include coordinates that equal one of the provided elements for the specified attribute or 
     * dimension. This method is not recommended for floating-point components.
     *
     * For monostate: no constraint is applied.
     *
     * @param col_index The index of the column 
     * @param selection The coordinate selection to apply to the columns.
     */
    SOMAIndexValueFilter& add_column_selection(int64_t col_index, SOMAColumnSelection<std::string> selection);

    /**
     * Returns `true` if at least one query condition is set.
     */
    bool is_initialized() const;

    /**
     * Returns a SOMA query condition that queries for the currently set coordinates.
     */
    SOMAValueFilter get_soma_query_condition() const;

    /**
     * Returns a TileDB query condition that queries for the currently set coordinates.
     * 
     * Throws an error if no query conditions are set.
     */
    inline QueryCondition get_query_condition() const {
        return get_soma_query_condition().query_condition();
    }

   private:
    /**Add a query condition to a coordinate. */
    SOMAIndexValueFilter& add_coordinate_query_condition(int64_t index, SOMAValueFilter&& qc);

    /** TileDB context required for creating query conditions. */
    std::shared_ptr<Context> ctx_;

    /** Names of the dimensions of the underlying TileDB Array. */
    std::vector<std::string> dim_names_;

    /** Query condition for coordinates. */
    std::vector<SOMAValueFilter> coord_qc_;
};

}  // namespace tiledbsoma

#endif
