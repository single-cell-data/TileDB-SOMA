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
        return QueryConditionExperimental::create<T>(ctx, elem_name, values, TILEDB_IN);
    }

    /**
     * Create a multi-condition QueryCondition using the `TILEDB_AND` predicate
     * to combine.
     *
     * For uninitialized query conditions:
     * * If neither query condition is initialized, this is left uninitialized.
     * * If only one object is initialized, this will be modified to be that condition.
     *
     * @param The other QueryCondition to add.
     * @returns A reference to this object for chaining methods.
     */
    SOMAQueryCondition& combine_with_and(const SOMAQueryCondition& other);

    /**
     * Create a multi-condition QueryCondition using the `TILEDB_AND` predicate
     * to combine.
     *
     * For uninitialized query conditions:
     * * If neither query condition is initialized, this is left uninitialized.
     * * If only one object is initialized, this will be modified to be that condition.
     *
     * @param The other QueryCondition to add.
     * @returns A reference to this object for chaining methods.
     */
    SOMAQueryCondition& combine_with_or(const SOMAQueryCondition& other);

    /**
     * Return if the query condition is initialized.
     */
    inline bool is_initialized() const {
        return qc_.has_value();
    }

    inline const QueryCondition& query_condition() const {
        return qc_.value();
    }

   private:
    /** Query condition for coordinates. */
    std::optional<QueryCondition> qc_;
};

class SOMACoordQueryCondition {
   public:
    SOMACoordQueryCondition() = delete;
    SOMACoordQueryCondition(const SOMAContext& ctx);

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
    SOMACoordQueryCondition& add_range(const std::string& elem_name, T start_value, T stop_value) {
        auto qc_range = QueryCondition::create<T>(*ctx_, elem_name, start_value, TILEDB_GE)
                            .combine(QueryCondition::create<T>(*ctx_, elem_name, stop_value, TILEDB_LE), TILEDB_AND);
        if (initialized_) {
            qc_ = qc_.combine(qc_range, TILEDB_AND);
        } else {
            qc_ = qc_range;
            initialized_ = true;
        }
        return *this;
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
    SOMACoordQueryCondition& add_points(const std::string& elem_name, const std::vector<T>& values) {
        auto qc_coords = QueryConditionExperimental::create<T>(*ctx_, elem_name, values, TILEDB_IN);
        if (initialized_) {
            qc_ = qc_.combine(qc_coords, TILEDB_AND);
        } else {
            qc_ = qc_coords;
            initialized_ = true;
        }
        return *this;
    }

    inline bool is_initialized() const {
        return initialized_;
    }

    inline const QueryCondition& query_condition() const {
        return qc_;
    }

   private:
    /** TileDB context required for creating query conditions. */
    std::shared_ptr<Context> ctx_;

    /** Query condition for coordinates. */
    QueryCondition qc_;

    /** Flag denoting it the query condition has been initialized yet. */
    bool initialized_ = false;
};

}  // namespace tiledbsoma

#endif
