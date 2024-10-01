#ifndef TILEDBSOMA_MULTIPOLYGON_H
#define TILEDBSOMA_MULTIPOLYGON_H

#include <vector>

#include "polygon.h"

namespace tiledbsoma
{
class MultiPolygon
{
public:
    MultiPolygon(std::vector<Polygon>&& polygons);
    ~MultiPolygon();

    std::vector<Polygon> polygons;
};
} // namespace tiledbsoma

#endif // TILEDBSOMA_MULTIPOLYGON_H
