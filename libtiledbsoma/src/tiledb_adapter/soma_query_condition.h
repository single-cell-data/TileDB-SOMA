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
#include <vector>

#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>
#include "../soma/soma_context.h"
#include "../utils/common.h"

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
     * Create a query condition from a range on an element.
     *
     * @tparam The type of the element the query condition applies to.
     * @param elem_name The name of the element the query condition applies to.
     * @param start_value The start of the range (inclusive).
     * @param stop_value The end of the range (inclusive).
     */
    template <typename T>
    static SOMAQueryCondition create_from_range(
        const Context& ctx, const std::string& elem_name, T start_value, T stop_value) {
        return SOMAQueryCondition(
            QueryCondition::create<T>(ctx, elem_name, start_value, TILEDB_GE)
                .combine(QueryCondition::create<T>(ctx, elem_name, stop_value, TILEDB_LE), TILEDB_AND));
    }

    /**
     * Create a query condition from a series of points on a range.
     *
     * @tparam The type of the element the query condition applies to.
     * @param elem_name The name of the element the query condition applies to.
     * @param values The values of the points
     */
    template <typename T>
    static SOMAQueryCondition create_from_points(
        const Context& ctx, const std::string& elem_name, const std::vector<T>& values) {
        if (values.empty()) {
            return SOMAQueryCondition();
        }
        return QueryConditionExperimental::create<T>(ctx, elem_name, values, TILEDB_IN);
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

class SOMACoordQueryCondition {
   public:
    SOMACoordQueryCondition() = delete;
    SOMACoordQueryCondition(const SOMAContext& ctx, const std::vector<std::string>& dim_names);

    SOMACoordQueryCondition(SOMACoordQueryCondition&& other) = default;
    SOMACoordQueryCondition(const SOMACoordQueryCondition& other) = default;
    ~SOMACoordQueryCondition() = default;

    /**
     * @brief Add a coordinate range to the query condition.
     *
     * The query condition will only select on elements of the TileDB array that include
     * coordinates within the selected range for the specified attribute or dimension.
     *
     * @tparam T The type of the column the query condition is applied to.
     * @param elem_name The name of the TileDB dimension or attribute the query condition is applied to.
     * @param start_value The first of the range.
     * @param stop_value The last value of the range.
     */
    template <typename T>
    SOMACoordQueryCondition& add_range(int64_t dim_index, T start_value, T stop_value) {
        return add_coordinate_query_condition(
            dim_index, SOMAQueryCondition::create_from_range<T>(*ctx_, dim_names_[dim_index], start_value, stop_value));
    }

    /**
     * @brief Add a list of coordinates to the query condition.
     *
     * The query condition will only select on elements of the TileDB array that include
     * coordinates that equal one of the provided elements for the specified attribute or 
     * dimension. This method is not recommended for floating-point components.
     *
     * Updates the query condition so only 
     * @tparam T The type of the column the query condition is applied to.
     * @param elem_name The name of the TileDB dimension or attribute the query condition is applied to.
     * @param values The values the coordinate will select on.
     *
     */
    template <typename T>
    SOMACoordQueryCondition& add_points(int64_t dim_index, const std::vector<T>& values) {
        return add_coordinate_query_condition(
            dim_index, SOMAQueryCondition::create_from_points<T>(*ctx_, dim_names_[dim_index], values));
    }

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
    QueryCondition get_query_condition() const {
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
