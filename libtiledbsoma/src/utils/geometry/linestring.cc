#include "linestring.h"

namespace tiledbsoma::geometry {
LineString::LineString(std::vector<BasePoint>&& points)
    : points(points) {
}

LineString::~LineString() {
}
}  // namespace tiledbsoma::geometry
