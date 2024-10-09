#include "write.h"

namespace tiledbsoma::geometry {

struct WKBSizeOperator {
    size_t binary_size(const BasePoint& point) {
        size_t size = 16;
        if (point.z.has_value())
            size += 8;
        if (point.m.has_value())
            size += 8;
        return size;
    }

    size_t operator()(const Point& point) {
        return 5 + binary_size(point);
    }

    size_t operator()(const LineString& linestring) {
        if (linestring.points.size() == 0)
            return 9;

        return 9 + linestring.points.size() *
                       binary_size(linestring.points.front());
    }

    size_t operator()(const Polygon& polygon) {
        size_t size = 1 + 4 + 4;

        if (polygon.exteriorRing.size() != 0)
            size += 4 + polygon.exteriorRing.size() *
                            binary_size(polygon.exteriorRing.front());
        else
            size += 4;

        for (auto& ring : polygon.interiorRings)
            if (ring.size() != 0)
                size += 4 + ring.size() * binary_size(ring.front());
            else
                size += 4;

        return size;
    }

    size_t operator()(const MultiPoint& multi_point) {
        size_t size = 9;

        for (auto& point : multi_point.points)
            size += this->operator()(point);

        return size;
    }

    size_t operator()(const MultiLineString& multi_linestring) {
        size_t size = 9;

        for (auto& linestring : multi_linestring.linestrings)
            size += this->operator()(linestring);

        return size;
    }

    size_t operator()(const MultiPolygon& multi_polygon) {
        size_t size = 9;

        for (auto& polygon : multi_polygon.polygons)
            size += this->operator()(polygon);

        return size;
    }

    size_t operator()(const GeometryCollection& collection) {
        size_t size = 9;

        for (auto& geometry : collection)
            size += std::visit(WKBSizeOperator{}, geometry);

        return size;
    }

    // size_t operator()(const auto& geometry) {
    //     throw std::runtime_error("Unsupported geometry type");
    // }
};

size_t wkb_size(const GenericGeometry& geometry) {
    return std::visit(WKBSizeOperator{}, geometry);
}

struct WKBWriteOperator {
    WKBWriteOperator(uint8_t* buffer, size_t& position, size_t size)
        : buffer(buffer)
        , position(position)
        , size(size) {
    }

    template <typename T>
    void write(const T& value) {
        assert(sizeof(T) + this->position <= this->size);

        memcpy(this->buffer + this->position, &value, sizeof(T));
        this->position += sizeof(T);
    }

    void wkb_write(const BasePoint& point) {
        write(point.x);
        write(point.y);
    }

    void operator()(const Point& point) {
        write((uint8_t)1);
        write(static_cast<uint32_t>(GeometryType::POINT));
        wkb_write(point);
    }

    void operator()(const LineString& linestring) {
        write((uint8_t)1);
        write(static_cast<uint32_t>(GeometryType::LINESTRING));
        write((uint32_t)linestring.points.size());

        for (auto& point : linestring.points) {
            wkb_write(point);
        }
    }

    void operator()(const Polygon& polygon) {
        write((uint8_t)1);
        write(static_cast<uint32_t>(GeometryType::POLYGON));
        write((uint32_t)(polygon.interiorRings.size() + 1));

        write((uint32_t)polygon.exteriorRing.size());
        for (auto& point : polygon.exteriorRing) {
            wkb_write(point);
        }

        for (auto& ring : polygon.interiorRings) {
            write((uint32_t)ring.size());
            for (auto& point : ring) {
                wkb_write(point);
            }
        }
    }

    void operator()(const MultiPoint& multi_point) {
        write((uint8_t)1);
        write(static_cast<uint32_t>(GeometryType::MULTIPOINT));
        write((uint32_t)multi_point.points.size());
        for (auto& point : multi_point.points)
            this->operator()(point);
    }

    void operator()(const MultiLineString& multi_linestring) {
        write((uint8_t)1);
        write(static_cast<uint32_t>(GeometryType::MULTILINESTRING));
        write((uint32_t)multi_linestring.linestrings.size());
        for (auto& linestring : multi_linestring.linestrings)
            this->operator()(linestring);
    }

    void operator()(const MultiPolygon& multi_polygon) {
        write((uint8_t)1);
        write(static_cast<uint32_t>(GeometryType::MULTIPOLYGON));
        write((uint32_t)multi_polygon.polygons.size());
        for (auto& polygon : multi_polygon.polygons)
            this->operator()(polygon);
    }

    void operator()(const GeometryCollection& collection) {
        write((uint8_t)1);
        write(static_cast<uint32_t>(GeometryType::GEOMETRYCOLLECTION));
        write((uint32_t)collection.size());
        for (auto& geometry : collection)
            std::visit(
                WKBWriteOperator{this->buffer, this->position, this->size},
                geometry);
    }

    void operator()(const auto&) {
        throw std::runtime_error("Unsupported geometry type");
    }

    uint8_t* buffer;
    size_t& position;
    size_t size;
};

void to_wkb(const GenericGeometry& geometry, uint8_t* buffer, size_t size) {
    size_t position = 0;

    std::visit(WKBWriteOperator{buffer, position, size}, geometry);

    assert(position == size);
}

BinaryBuffer to_wkb(const GenericGeometry& geometry) {
    BinaryBuffer buffer(wkb_size(geometry));
    to_wkb(geometry, buffer.data(), buffer.size());
    return buffer;
}

}  // namespace tiledbsoma::geometry