#include "point.h"

namespace tiledbsoma::geometry {
Point::Point()
    : BasePoint(0, 0){};
Point::Point(
    double_t x,
    double_t y,
    std::optional<double_t> z,
    std::optional<double_t> m)
    : BasePoint(x, y, z, m) {
}
Point::Point(BasePoint&& point)
    : BasePoint(std::move(point)){};

Point::~Point() {
}
}  // namespace tiledbsoma::geometry
