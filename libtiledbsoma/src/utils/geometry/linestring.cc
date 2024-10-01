#include "linestring.h"

namespace tiledbsoma
{
    LineString::LineString(std::vector<BasePoint>&& points) : points(points) {}

    LineString::~LineString() {}
} // namespace tiledbsoma
