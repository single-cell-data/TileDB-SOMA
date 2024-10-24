#ifndef TILEDBSOMA_GEOMETRY_H
#define TILEDBSOMA_GEOMETRY_H

#include <variant>

#include "linestring.h"
#include "multilinestring.h"
#include "multipoint.h"
#include "multipolygon.h"
#include "point.h"
#include "polygon.h"

namespace tiledbsoma::geometry {

enum GeometryType : uint32_t {
    POINT = 1,
    LINESTRING = 2,
    POLYGON = 3,
    MULTIPOINT = 4,
    MULTILINESTRING = 5,
    MULTIPOLYGON = 6,
    GEOMETRYCOLLECTION = 7
};

using BinaryBuffer = std::vector<std::byte>;

struct GeometryCollection;
using GenericGeometry = std::variant<
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection>;

struct GeometryCollection : public std::vector<GenericGeometry> {
    using std::vector<GenericGeometry>::vector;
};
}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_GEOMETRY_H