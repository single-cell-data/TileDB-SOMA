#include "multipoint.h"

namespace tiledbsoma::geometry {
MultiPoint::MultiPoint(std::vector<Point>&& points)
    : points(points) {
}

MultiPoint::~MultiPoint() {
}
}  // namespace tiledbsoma::geometry
