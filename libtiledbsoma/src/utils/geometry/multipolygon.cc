#include "multipolygon.h"

namespace tiledbsoma
{
    MultiPolygon::MultiPolygon(std::vector<Polygon>&& polygons) : polygons(polygons) {}

    MultiPolygon::~MultiPolygon() {}
} // namespace tiledbsoma
