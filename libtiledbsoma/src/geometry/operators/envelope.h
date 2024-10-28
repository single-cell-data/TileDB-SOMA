#ifndef TILEDBSOMA_GEOMETRY_ENVELOPE_H
#define TILEDBSOMA_GEOMETRY_ENVELOPE_H

#include <array>

#include "../geometry.h"

namespace tiledbsoma::geometry {

struct Envelope {
    Envelope();

    std::array<std::pair<double_t, double_t>, 4> range;
};

struct EnvelopeOperator {
    EnvelopeOperator(Envelope& envelope);

    void base_envelope(const BasePoint& point);
    void operator()(const Point& point);
    void operator()(const LineString& linestring);
    void operator()(const Polygon& polygon);
    void operator()(const MultiPoint& multi_point);
    void operator()(const MultiLineString& multi_linestring);
    void operator()(const MultiPolygon& multi_polygon);
    void operator()(const GeometryCollection& collection);

    Envelope& envelope;
};

Envelope envelope(const GenericGeometry& geometry);

}  // namespace tiledbsoma::geometry
#endif