#include "multipoint.h"

namespace tiledbsoma
{
    MultiPoint::MultiPoint(std::vector<Point>&& points) : points(points) {}

    MultiPoint::~MultiPoint() {}
} // namespace tiledbsoma
