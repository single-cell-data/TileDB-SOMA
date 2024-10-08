#ifndef TILEDBSOMA_GEOMETRY_H
#define TILEDBSOMA_GEOMETRY_H

#include <variant>

#include "point.h"
#include "linestring.h"
#include "polygon.h"
#include "multipoint.h"
#include "multilinestring.h"
#include "multipolygon.h"

namespace tiledbsoma::geometry {

enum class GeometryType {
    POINT = 1,
    LINESTRING = 2,
    POLYGON = 3,
    MULTIPOINT = 4,
    MULTILINESTRING = 5,
    MULTIPOLYGON = 6,
    GEOMETRYCOLLECTION = 7
};

using BinaryBuffer = std::vector<uint8_t>;

struct GeometryCollection;
using GenericGeometry = std::variant<Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection>;

struct GeometryCollection : public std::vector<GenericGeometry>
{
    using std::vector<GenericGeometry>::vector;
};
}

#endif  // TILEDBSOMA_GEOMETRY_H