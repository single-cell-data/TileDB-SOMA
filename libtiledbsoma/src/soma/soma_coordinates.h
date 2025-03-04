/** * @file   soma_coordinates.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines classes, structs, and helpers for managing coordinate
 *   spaces and coordinate space transformations.
 */

#ifndef SOMA_METADATA
#define SOMA_METADATA

#include <optional>
#include <string>
#include <vector>
#include "../utils/common.h"

namespace tiledbsoma {

struct SOMAAxis {
    std::string name{"x"};
    std::optional<std::string> unit{std::nullopt};

    inline friend bool operator==(const SOMAAxis& lhs, const SOMAAxis& rhs) {
        return lhs.name == rhs.name && lhs.unit == rhs.unit;
    }

    inline friend bool operator!=(const SOMAAxis& lhs, const SOMAAxis& rhs) {
        return !(lhs == rhs);
    }
};

class SOMACoordinateSpace {
   public:
    static SOMACoordinateSpace from_metadata(
        tiledb_datatype_t value_type, uint32_t value_num, const void* value);

    static SOMACoordinateSpace from_string(std::string_view metadata);

    SOMACoordinateSpace();

    SOMACoordinateSpace(const std::vector<SOMAAxis>& axes);

    SOMACoordinateSpace(const std::vector<std::string>& axis_names);

    SOMACoordinateSpace(
        const std::vector<std::string>& axis_names,
        const std::vector<std::optional<std::string>>& axis_units);

    inline friend bool operator==(
        const SOMACoordinateSpace& lhs, const SOMACoordinateSpace& rhs) {
        return lhs.axes_ == rhs.axes_;
    }

    inline friend bool operator!=(
        const SOMACoordinateSpace& lhs, const SOMACoordinateSpace& rhs) {
        return !(lhs == rhs);
    }

    std::string to_string() const;

    inline size_t size() const {
        return axes_.size();
    }

    inline const SOMAAxis& axis(size_t index) const {
        return axes_[index];
    }

   private:
    std::vector<SOMAAxis> axes_;
};

}  // namespace tiledbsoma

#endif
