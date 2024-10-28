#include "envelope.h"

#include <limits>

namespace tiledbsoma::geometry {
Envelope::Envelope() {
    this->range = {
        std::make_pair(
            std::numeric_limits<double_t>::max(),
            -std::numeric_limits<double_t>::max()),
        std::make_pair(
            std::numeric_limits<double_t>::max(),
            -std::numeric_limits<double_t>::max()),
        std::make_pair(
            std::numeric_limits<double_t>::max(),
            -std::numeric_limits<double_t>::max()),
        std::make_pair(
            std::numeric_limits<double_t>::max(),
            -std::numeric_limits<double_t>::max())};
}

EnvelopeOperator::EnvelopeOperator(Envelope& envelope)
    : envelope(envelope) {
}

void EnvelopeOperator::base_envelope(const BasePoint& point) {
    this->envelope.range[0].first = std::min(
        this->envelope.range[0].first, point.x);
    this->envelope.range[0].second = std::max(
        this->envelope.range[0].second, point.x);

    this->envelope.range[1].first = std::min(
        this->envelope.range[1].first, point.y);
    this->envelope.range[1].second = std::max(
        this->envelope.range[1].second, point.y);

    if (point.z.has_value()) {
        this->envelope.range[2].first = std::min(
            this->envelope.range[2].first, point.z.value());
        this->envelope.range[2].second = std::max(
            this->envelope.range[2].second, point.z.value());
    }

    if (point.m.has_value()) {
        this->envelope.range[3].first = std::min(
            this->envelope.range[3].first, point.m.value());
        this->envelope.range[3].second = std::max(
            this->envelope.range[3].second, point.m.value());
    }
}

void EnvelopeOperator::operator()(const Point& point) {
    this->base_envelope(point);
}

void EnvelopeOperator::operator()(const LineString& linestring) {
    for (auto& point : linestring.points) {
        this->base_envelope(point);
    }
}

void EnvelopeOperator::operator()(const Polygon& polygon) {
    for (auto& point : polygon.exteriorRing) {
        base_envelope(point);
    }
}

void EnvelopeOperator::operator()(const MultiPoint& multi_point) {
    for (auto& point : multi_point.points) {
        this->operator()(point);
    }
}

void EnvelopeOperator::operator()(const MultiLineString& multi_linestring) {
    for (auto& linestring : multi_linestring.linestrings) {
        this->operator()(linestring);
    }
}

void EnvelopeOperator::operator()(const MultiPolygon& multi_polygon) {
    for (auto& polygon : multi_polygon.polygons) {
        this->operator()(polygon);
    }
}

void EnvelopeOperator::operator()(const GeometryCollection& collection) {
    for (auto& geometry : collection) {
        std::visit(EnvelopeOperator{this->envelope}, geometry);
    }
}

Envelope envelope(const GenericGeometry& geometry) {
    Envelope env;

    std::visit(EnvelopeOperator{env}, geometry);

    return env;
}
}  // namespace tiledbsoma::geometry
