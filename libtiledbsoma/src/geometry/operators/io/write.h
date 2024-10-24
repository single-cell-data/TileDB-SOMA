#ifndef TILEDBSOMA_GEOMETRY_WRITE_H
#define TILEDBSOMA_GEOMETRY_WRITE_H

#include <cassert>
#include <cstring>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "../../geometry.h"

namespace tiledbsoma::geometry {

const size_t WKB_BYTE_ORDER_SIZE = 1;
const size_t WKB_GEOEMTRY_TYPE_SIZE = 4;
const size_t WKB_ELEMENT_COUNT_SIZE = 4;

struct WKBSizeOperator {
    size_t binary_size(const BasePoint& point);
    size_t operator()(const Point& point);
    size_t operator()(const LineString& linestring);
    size_t operator()(const Polygon& polygon);
    size_t operator()(const MultiPoint& multi_point);
    size_t operator()(const MultiLineString& multi_linestring);
    size_t operator()(const MultiPolygon& multi_polygon);
    size_t operator()(const GeometryCollection& collection);
};

struct WKBWriteOperator {
    WKBWriteOperator(std::byte* buffer, size_t& position, size_t size);

    template <typename T>
    void write(const T& value) {
        assert(sizeof(T) + this->position <= this->size);

        memcpy(this->buffer + this->position, &value, sizeof(T));
        this->position += sizeof(T);
    }

    void wkb_write(const BasePoint& point);
    void operator()(const Point& point);
    void operator()(const LineString& linestring);
    void operator()(const Polygon& polygon);
    void operator()(const MultiPoint& multi_point);
    void operator()(const MultiLineString& multi_linestring);
    void operator()(const MultiPolygon& multi_polygon);
    void operator()(const GeometryCollection& collection);

    std::byte* buffer;
    size_t& position;
    size_t size;
};

size_t wkb_size(const GenericGeometry& geometry);

void to_wkb(const GenericGeometry& geometry, uint8_t* buffer, size_t size);

BinaryBuffer to_wkb(const GenericGeometry& geometry);

}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_GEOMETRY_WRITE_H