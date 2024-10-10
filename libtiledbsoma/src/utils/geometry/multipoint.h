#ifndef TILEDBSOMA_MULTIPOINT_H
#define TILEDBSOMA_MULTIPOINT_H

#include <vector>

#include "point.h"

namespace tiledbsoma::geometry {
class MultiPoint {
   public:
    MultiPoint(std::vector<Point>&& points = std::vector<Point>());
    ~MultiPoint();

    std::vector<Point> points;
};

}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_MULTIPOINT_H