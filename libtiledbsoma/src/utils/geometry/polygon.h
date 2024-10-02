#ifndef TILEDBSOMA_POLYGON_H
#define TILEDBSOMA_POLYGON_H

#include <vector>

#include "base.h"
#include "point.h"

namespace tiledbsoma
{
    class Polygon
    {
    public:
        Polygon(std::vector<BasePoint>&& exteriorRing = std::vector<BasePoint>(), std::vector<std::vector<BasePoint>>&& interiorRings = std::vector<std::vector<BasePoint>>());
        ~Polygon();

        std::vector<BasePoint> exteriorRing;
        std::vector<std::vector<BasePoint>> interiorRings;
    };
    
} // namespace tiledbsoma

#endif // TILEDBSOMA_POLYGON_H