#ifndef TILEDBSOMA_LINESTRING_H
#define TILEDBSOMA_LINESTRING_H

#include <vector>

#include "base.h"

namespace tiledbsoma::geometry {

class LineString {
   public:
    LineString(std::vector<BasePoint>&& points = std::vector<BasePoint>());
    ~LineString();

    std::vector<BasePoint> points;
};
}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_LINESTRING_H