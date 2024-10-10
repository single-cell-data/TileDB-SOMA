#include "polygon.h"

namespace tiledbsoma::geometry {
Polygon::Polygon(
    std::vector<BasePoint>&& exteriorRing,
    std::vector<std::vector<BasePoint>>&& interiorRings)
    : exteriorRing(exteriorRing)
    , interiorRings(interiorRings) {
}

Polygon::~Polygon() {
}
}  // namespace tiledbsoma::geometry
