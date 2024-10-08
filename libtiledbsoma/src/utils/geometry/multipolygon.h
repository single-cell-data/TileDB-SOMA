#ifndef TILEDBSOMA_MULTIPOLYGON_H
#define TILEDBSOMA_MULTIPOLYGON_H

#include <vector>

#include "polygon.h"

namespace tiledbsoma::geometry {
class MultiPolygon {
   public:
    MultiPolygon(std::vector<Polygon>&& polygons = std::vector<Polygon>());
    ~MultiPolygon();

    std::vector<Polygon> polygons;
};
}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_MULTIPOLYGON_H
