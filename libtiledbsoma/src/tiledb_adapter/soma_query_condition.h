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
 *
 * WARNING: Do not include in public facing header: this directly imports the logger
 * for fmt.
 */

#ifndef TILEDBSOMA_TILEDB_ADAPTER_H
#define TILEDBSOMA_TILEDB_ADAPTER_H

#include <memory>
#include <optional>
#include <span>
#include <vector>

#include <tiledb/tiledb.h>
#include <tiledb/tiledb>

#include "../common/soma_column_selection.h"
#include "../soma/soma_context.h"
#include "../utils/common.h"
#include "../utils/logger.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAQueryCondition {
   public:
    SOMAQueryCondition() = default;
    SOMAQueryCondition(SOMAQueryCondition&& other) = default;
    SOMAQueryCondition(const SOMAQueryCondition& other) = default;
    ~SOMAQueryCondition() = default;

    SOMAQueryCondition& operator=(const SOMAQueryCondition& other) = default;
    SOMAQueryCondition& operator=(SOMAQueryCondition&& other) = default;

    /**
     * Create a SOMAQueryCondition from an initialized TileDB query condition.
     *
     * Important: The TileDB QueryCondtion must be initialized.
     *
     * @param qc The initialized TileDB query condition.
     */
    SOMAQueryCondition(const QueryCondition& qc);

    /**
     * Create a SOMAQueryCondition from an initialized TileDB query condition.
     *
     * Important: The TileDB QueryCondtion must be initialized.
     *
     * @param qc The initialized TileDB query condition.
     */
    SOMAQueryCondition(QueryCondition&& qc);

    /**
     * Create a query condition from a slice.
     *
     * @tparam The type of the element the query condition applies to.
     * @param elem_name The name of the element the query condition applies to.
     * @param start_value The start of the range (inclusive).
     * @param stop_value The end of the range (inclusive).
     */
    template <typename T, typename = std::enable_if<!std::is_same_v<T, std::string>>>
    static SOMAQueryCondition create_from_range(
        const Context& ctx, const std::string& elem_name, T start_value, T stop_value) {
        if (stop_value < start_value) {
            throw std::invalid_argument(
                fmt::format(
                    "Cannot set range [{}, {}] on column '{}'. Invalid range: the final value must be greater "
                    "than or equal to the starting value.",
                    start_value,
                    stop_value,
                    elem_name));
        }
        return SOMAQueryCondition(
            QueryCondition::create<T>(ctx, elem_name, start_value, TILEDB_GE)
                .combine(QueryCondition::create<T>(ctx, elem_name, stop_value, TILEDB_LE), TILEDB_AND));
    }

    /**
     * Create a query condition from a slice.
     *
     * @param elem_name The name of the element the query condition applies to.
     * @param start_value The start of the range (inclusive).
     * @param stop_value The end of the range (inclusive).
     */
    static SOMAQueryCondition create_from_range(
        const Context& ctx,
        const std::string& elem_name,
        const std::string& start_value,
        const std::string& stop_value);

    /**
     * Create a query condition from a series of points on a range.
     *
     * @tparam The type of the element the query condition applies to.
     * @param elem_name The name of the element the query condition applies to.
     * @param values The values of the points
     */
    template <typename T, typename = std::enable_if<!std::is_same_v<T, std::string>>>
    static SOMAQueryCondition create_from_points(
        const Context& ctx, const std::string& elem_name, std::span<T> values) {
        if (values.empty()) {
            throw std::invalid_argument(
                fmt::format("Cannot set coordinates on column '{}'. No coordinates provided.", elem_name));
        }
        // Using C API because C++ API only supports std::vector, not std::span.
        std::vector<uint64_t> offsets(values.size());
        for (size_t index = 0; index < values.size(); ++index) {
            offsets[index] = index * sizeof(T);
        }
        tiledb_query_condition_t* qc;
        ctx.handle_error(tiledb_query_condition_alloc_set_membership(
            ctx.ptr().get(),
            elem_name.c_str(),
            values.data(),
            values.size() * sizeof(T),
            offsets.data(),
            offsets.size() * sizeof(uint64_t),
            TILEDB_IN,
            &qc));
        return QueryCondition(ctx, qc);
    }

    static SOMAQueryCondition create_from_points(
        const Context& ctx, const std::string& elem_name, std::span<std::string> values);
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

class SOMACoordQueryCondition {
   public:
    SOMACoordQueryCondition() = delete;
    SOMACoordQueryCondition(const SOMAContext& ctx, const std::vector<std::string>& dim_names);

    SOMACoordQueryCondition(SOMACoordQueryCondition&& other) = default;
    SOMACoordQueryCondition(const SOMACoordQueryCondition& other) = default;
    ~SOMACoordQueryCondition() = default;

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
    SOMACoordQueryCondition& add_column_selection(
        int64_t col_index, SOMAColumnSelection<T> selection, const std::pair<T, T>& domain) {
        std::visit(
            [&](auto&& val) {
                using S = std::decay_t<decltype(val)>;
                if constexpr (std::is_same_v<S, SOMASliceSelection<T>>) {
                    if (!val.has_overlap(domain)) {
                        throw std::out_of_range(
                            fmt::format(
                                "Non-overlapping slice [{}, {}] on column '{}'. Slice must overlap the current "
                                "column domain [{}, {}].",
                                val.start,
                                val.stop,
                                dim_names_[col_index],
                                domain.first,
                                domain.second));
                    }
                    add_coordinate_query_condition(
                        col_index,
                        SOMAQueryCondition::create_from_range<T>(*ctx_, dim_names_[col_index], val.start, val.stop));

                } else if constexpr (std::is_same_v<S, SOMAPointSelection<T>>) {
                    if (!val.is_subset(domain)) {
                        throw std::out_of_range(
                            fmt::format(
                                "Out-of-bounds coordinates found on column '{}'. Coordinates must be "
                                "inside column daomain [{}, {}].",
                                dim_names_[col_index],
                                domain.first,
                                domain.second));
                    }
                    add_coordinate_query_condition(
                        col_index, SOMAQueryCondition::create_from_points<T>(*ctx_, dim_names_[col_index], val.points));
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
    SOMACoordQueryCondition& add_column_selection(int64_t col_index, SOMAColumnSelection<std::string> selection);

    /**
     * Returns `true` if at least one query condition is set.
     */
    bool is_initialized() const;

    /**
     * Returns a SOMA query condition that queries for the currently set coordinates.
     */
    SOMAQueryCondition get_soma_query_condition() const;

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
    SOMACoordQueryCondition& add_coordinate_query_condition(int64_t index, SOMAQueryCondition&& qc);

    /** TileDB context required for creating query conditions. */
    std::shared_ptr<Context> ctx_;

    /** Names of the dimensions of the underlying TileDB Array. */
    std::vector<std::string> dim_names_;

    /** Query condition for coordinates. */
    std::vector<SOMAQueryCondition> coord_qc_;
};

}  // namespace tiledbsoma

#endif
