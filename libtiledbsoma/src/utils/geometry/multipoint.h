#ifndef TILEDBSOMA_MULTIPOINT_H
#define TILEDBSOMA_MULTIPOINT_H

#include <vector>

#include "point.h"

namespace tiledbsoma
{
    class MultiPoint
    {
    public:
        MultiPoint(std::vector<Point>&& points);
        ~MultiPoint();

        std::vector<Point> points;
    };
    
} // namespace tiledbsoma

#endif // TILEDBSOMA_MULTIPOINT_H