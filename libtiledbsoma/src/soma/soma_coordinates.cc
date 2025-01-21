/**
 *   This file defines classes, structs, and helpers for managing coordinate
 *   spaces and coordinate space transformations.
 * @file   soma_coordinates.cc
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

#include <tiledb/tiledb>
#include <unordered_set>
#include "../utils/logger.h"
#include "nlohmann/json.hpp"

#include "soma_coordinates.h"

using json = nlohmann::json;

namespace nlohmann {
template <>
struct adl_serializer<tiledbsoma::SOMAAxis> {
    static void to_json(json& j, const tiledbsoma::SOMAAxis& axis) {
        // void to_json(json& j, tiledbsoma::SOMAAxis& axis) {
        if (axis.unit.has_value()) {
            j = json{{"name", axis.name}, {"unit", axis.unit.value()}};
        } else {
            j = json{{"name", axis.name}, {"unit", nullptr}};
        }
    }

    static void from_json(const json& j, tiledbsoma::SOMAAxis& axis) {
        j.at("name").get_to(axis.name);
        auto unit_json = j.at("unit");
        if (unit_json.is_null()) {
            axis.unit = std::nullopt;
        } else {
            unit_json.get_to(axis.unit);
        }
    }
};
}  // namespace nlohmann

namespace tiledbsoma {

SOMACoordinateSpace::SOMACoordinateSpace()
    : axes_{{"x", std::nullopt}, {"y", std::nullopt}} {
}

SOMACoordinateSpace::SOMACoordinateSpace(const std::vector<SOMAAxis>& axes)
    : axes_{axes} {
    if (axes_.size() == 0) {
        throw TileDBSOMAError("Coordinate space must have at least one axis.");
    }
    std::unordered_set<std::string> axis_names;
    for (const auto& axis : axes_) {
        axis_names.emplace(axis.name);
    }
    if (axes_.size() != axis_names.size()) {
        throw TileDBSOMAError(
            "The name for coordinate space axes must be unique.");
    }
}

SOMACoordinateSpace::SOMACoordinateSpace(
    const std::vector<std::string>& axis_names) {
    if (axis_names.size() == 0) {
        throw TileDBSOMAError("Coordinate space must have at least one axis.");
    }
    std::unordered_set<std::string> unique_axis_names(
        axis_names.begin(), axis_names.end());
    if (axis_names.size() != unique_axis_names.size()) {
        throw TileDBSOMAError(
            "The name for coordinate space axes must be unique.");
    }
    axes_.reserve(axis_names.size());
    for (const auto& name : axis_names) {
        axes_.push_back({name, std::nullopt});
    }
}

SOMACoordinateSpace::SOMACoordinateSpace(
    const std::vector<std::string>& axis_names,
    const std::vector<std::optional<std::string>>& axis_units) {
    if (axis_names.size() != axis_units.size()) {
        throw TileDBSOMAError(
            "[SOMACoordinateSpace]: Axis names and axis units size mismatch. ");
    }
    auto num_axes = axis_names.size();
    if (num_axes == 0) {
        throw TileDBSOMAError("Coordinate space must have at least one axis.");
    }
    std::unordered_set<std::string> unique_axis_names(
        axis_names.begin(), axis_names.end());
    if (axis_names.size() != unique_axis_names.size()) {
        throw TileDBSOMAError(
            "The name for coordinate space axes must be unique.");
    }
    axes_.reserve(num_axes);
    for (size_t index{0}; index < num_axes; ++index) {
        axes_.push_back({axis_names[index], axis_units[index]});
    }
}

SOMACoordinateSpace SOMACoordinateSpace::from_metadata(
    tiledb_datatype_t value_type, uint32_t value_num, const void* value) {
    if (value_type != TILEDB_STRING_UTF8 && value_type != TILEDB_STRING_ASCII) {
        throw TileDBSOMAError(fmt::format(
            "[SOMACoordinateSpace]: Unexpected datatype for coordinate space "
            "metadata. Expected {} or {}; got {}",
            tiledb::impl::type_to_str(TILEDB_STRING_UTF8),
            tiledb::impl::type_to_str(TILEDB_STRING_ASCII),
            tiledb::impl::type_to_str(value_type)));
    }
    if (value == nullptr) {
        throw TileDBSOMAError(
            "[SOMACoordinateSpace]: Missing value for coordinate space "
            "metadata.");
    }
    std::string value_str(static_cast<const char*>(value), value_num);
    auto value_json = json::parse(value_str);
    auto axes = value_json.template get<std::vector<SOMAAxis>>();

    return SOMACoordinateSpace(axes);
}

std::string SOMACoordinateSpace::to_string() const {
    json serializer(axes_);
    return serializer.dump(-1, ' ', true);
}
}  // namespace tiledbsoma
