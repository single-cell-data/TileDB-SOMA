/**
 * @file   soma_coordinate_selection.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This defines the coordinate seletion API.
 *
 */

#ifndef SOMA_COORDINATE_SELECTION_H
#define SOMA_COORDINATE_SELECTION_H

#include <math.h>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <variant>

#include "../utils/common.h"

namespace tiledbsoma {

template <typename T>
struct SOMASliceSelection {
   public:
    /**
     * Creates a closed slice (includes both end points).
     *
     * @param slice_start The first value (inclusive) of the slice.
     * @param slice_stop The last value (inclusive) of the slice.
     */
    SOMASliceSelection(T slice_start, T slice_stop)
        : start{slice_start}
        , stop{slice_stop} {
        // Using sstream because we don't want to include fmt directly in external header.
        std::stringstream ss;
        ss << "Invalid slice [ " << start << ", " << stop
           << "]. The lower bound must be less than or equal to the upper bound.";
        if (stop < start) {
            throw std::invalid_argument(ss.str());
        }
    }

    SOMASliceSelection(const std::pair<T, T>& slice)
        : SOMASliceSelection(slice.first, slice.second) {};

    /** Returns if the slice overlaps a requested interval. 
     *
     * The check assumes both the slice objet and the interval are closed intervals
     * that include the end points. The provided interval must be a valid range with 
     * `interval.first <= interval.second`.
     *
     * @param interval The interval to check overlap against.
     */
    bool has_overlap(std::pair<T, T> interval) {
        return (stop >= interval.first && start <= interval.second);
    }

    T start;
    T stop;
};

template <typename T>
struct SOMAPointSelection {
   public:
    SOMAPointSelection(std::span<T> point_data)
        : points{point_data} {
    }

    /** Returns if the points are strictly contained within the requested interval. */
    bool is_subset(const std::pair<T, T>& interval) {
        /** Can switch to the following for C++23
        return std::any_of(
            points.cbegin(), points.cend(), [&](auto val) { return val < interval.first || val > interval.second; });
        */
        for (const auto& val : points) {
            if (val < interval.first || val > interval.second) {
                return false;
            }
        }
        return true;
    }

    /** The points to select for. */
    std::span<T> points;
};

template <typename T>
using SOMAColumnSelection = std::variant<std::monostate, SOMASliceSelection<T>, SOMAPointSelection<T>>;

}  // namespace tiledbsoma

#endif
