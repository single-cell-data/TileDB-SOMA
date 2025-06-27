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
#include <tiledb/tiledb>
#include <tiledb/tiledb_experimental>
#include "../soma/soma_context.h"
#include "../utils/common.h"

namespace tiledbsoma {
using namespace tiledb;

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
    void add_range(const std::string& elem_name, T start_value, T stop_value) {
        auto qc_range = QueryCondition::create<T>(*ctx_, elem_name, start_value, TILEDB_GE)
                            .combine(QueryCondition::create<T>(*ctx_, elem_name, stop_value, TILEDB_LE), TILEDB_AND);
        if (initialized_) {
            qc_ = qc_.combine(qc_range, TILEDB_AND);
        } else {
            qc_ = qc_range;
            initialized_ = true;
        }
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
    void add_points(const std::string& elem_name, const std::vector<T>& values) {
        auto qc_coords = QueryConditionExperimental::create<T>(*ctx_, elem_name, values, TILEDB_IN);
        if (initialized_) {
            qc_ = qc_.combine(qc_coords, TILEDB_AND);
        } else {
            qc_ = qc_coords;
            initialized_ = true;
        }
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
