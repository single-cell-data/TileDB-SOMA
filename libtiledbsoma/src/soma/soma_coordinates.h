/**
 * @file   soma_coordinates.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2025 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
