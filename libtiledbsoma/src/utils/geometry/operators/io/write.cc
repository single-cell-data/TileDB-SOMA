#include "write.h"

namespace tiledbsoma::geometry
{
    size_t wkb_size(const GenericGeometry& geometry) {
        auto binary_size = [](const BasePoint& point) -> size_t {
            size_t size = 16;
            if (point.z.has_value()) size += 8;
            if (point.m.has_value()) size += 8;
            return size;
        };

        Operator compute_size = {
            [&](auto&, const Point& point) -> size_t {
                return 5 + binary_size(point);
            },
            [&](auto&, const LineString& linestring) -> size_t {
                if (linestring.points.size() == 0) return 9;

                return 9 + linestring.points.size() * binary_size(linestring.points.front());
            },
            [&](auto&, const Polygon& polygon) -> size_t {
                size_t size =  1 + 4 + 4;

                if (polygon.exteriorRing.size() != 0)
                    size += 4 + polygon.exteriorRing.size() * binary_size(polygon.exteriorRing.front());
                else
                    size += 4;

                for (auto& ring : polygon.interiorRings)
                    if (ring.size() != 0) 
                        size += 4 + ring.size() * binary_size(ring.front());
                    else 
                        size += 4;

                return size;
            },
            [&](auto& self, const MultiPoint& multi_point) -> size_t {
                size_t size = 9;

                for (auto& point : multi_point.points)
                    size += std::visit(self, GenericGeometry(point));

                return size;
            },
            [&](auto& self, const MultiLineString& multi_linestring) -> size_t {
                size_t size = 9;

                for (auto& linestring : multi_linestring.linestrings)
                    size += std::visit(self, GenericGeometry(linestring));

                return size;
            },
            [&](auto& self, const MultiPolygon& multi_polygon) -> size_t {
                size_t size = 9;

                for (auto& polygon : multi_polygon.polygons)
                    size += std::visit(self, GenericGeometry(polygon));

                return size;
            },
            [&](auto& self, const GeometryCollection& collection) -> size_t {
                size_t size = 9;

                for (auto& geometry : collection)
                    size += std::visit(self, geometry);

                return size;
            },
            [&](auto, auto) -> size_t {
                throw std::runtime_error("Unsupported geometry type");
            }
        };

        return std::visit(Recurse{compute_size}, geometry);
    }

    void to_wkb(const GenericGeometry& geometry, uint8_t* buffer, size_t size) {
        size_t position = 0;

        auto write = [&buffer, &position, &size]<typename T>(const T& value) -> void {
            assert(sizeof(T) + position <= size);

            memcpy(buffer + position, &value, sizeof(T));
            position += sizeof(T);
        };

        Operator serialize {
            [&](auto&, const Point& point) -> void {
                write((uint8_t)1);
                write(static_cast<uint32_t>(GeometryType::POINT));
                write(point.x);
                write(point.y);
            },
            [&](auto&, const LineString& linestring) -> void {
                write((uint8_t)1);
                write(static_cast<uint32_t>(GeometryType::LINESTRING));
                write((uint32_t)linestring.points.size());

                for (auto& point : linestring.points) {
                    write(point.x);
                    write(point.y);
                }
            },
            [&](auto&, const Polygon& polygon) -> void {
                write((uint8_t)1);
                write(static_cast<uint32_t>(GeometryType::POLYGON));
                write((uint32_t)(polygon.interiorRings.size() + 1));

                write((uint32_t)polygon.exteriorRing.size());
                for (auto& point : polygon.exteriorRing) {
                    write(point.x);
                    write(point.y);
                }

                for (auto& ring : polygon.interiorRings) {
                    write((uint32_t)ring.size());
                    for (auto& point : ring) {
                        write(point.x);
                        write(point.y);
                    }
                }
            },
            [&](auto& self, const MultiPoint& multi_point) -> void {
                write((uint8_t)1);
                write(static_cast<uint32_t>(GeometryType::MULTIPOINT));
                write((uint32_t)multi_point.points.size());
                for (auto& point : multi_point.points)
                    std::visit(self, GenericGeometry(point));
            },
            [&](auto& self, const MultiLineString& multi_linestring) -> void {
                write((uint8_t)1);
                write(static_cast<uint32_t>(GeometryType::MULTILINESTRING));
                write((uint32_t)multi_linestring.linestrings.size());
                for (auto& linestring : multi_linestring.linestrings)
                    std::visit(self, GenericGeometry(linestring));
            },
            [&](auto& self, const MultiPolygon& multi_polygon) -> void {
                write((uint8_t)1);
                write(static_cast<uint32_t>(GeometryType::MULTIPOLYGON));
                write((uint32_t)multi_polygon.polygons.size());
                for (auto& polygon : multi_polygon.polygons)
                    std::visit(self, GenericGeometry(polygon));
            },
            [&](auto& self, const GeometryCollection& collection) -> void {
                write((uint8_t)1);
                write(static_cast<uint32_t>(GeometryType::GEOMETRYCOLLECTION));
                write((uint32_t)collection.size());
                for (auto& geometry : collection)
                    std::visit(self, geometry);
            },
            [](auto, auto) -> void {
                throw std::runtime_error("Unsupported geometry type");
            }
        };

        std::visit(Recurse{serialize}, geometry);
    }

    BinaryBuffer to_wkb(const GenericGeometry& geometry) {
        BinaryBuffer buffer(wkb_size(geometry));
        to_wkb(geometry, buffer.data(), buffer.size());
        return buffer;
    }

    
} // namespace tiledbsoma::geometry