#include "multipolygon.h"

namespace tiledbsoma::geometry {
MultiPolygon::MultiPolygon(std::vector<Polygon>&& polygons)
    : polygons(polygons) {
}

MultiPolygon::~MultiPolygon() {
}
}  // namespace tiledbsoma::geometry
