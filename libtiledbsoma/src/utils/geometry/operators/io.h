#ifndef TILEDBSOMA_GEOMETRY_IO_H
#define TILEDBSOMA_GEOMETRY_IO_H

#include <vector>
#include <string>
#include <cstring>
#include <variant>
#include <stdexcept>
#include <cassert>

#include "../geometry.h"

namespace tiledbsoma
{
    template <typename S> struct Reader {};
    
    template <>
    struct Reader<BinaryBuffer>
    {
        Reader(BinaryBuffer& buffer) : buffer(buffer), position(0) {}

        template <typename T> T read() {
            assert(this->position + sizeof(T) <= this->buffer.size());

            T value = *(T*)(&this->buffer[this->position]);
            this->position += sizeof(T);
            return value;
        }

        template <typename T> T peek(size_t offset = 0) const {
            assert(this->position + offset + sizeof(T) <= this->buffer.size());

            T value = *(T*)(&this->buffer[this->position + offset]);
            return value;
        }

        std::vector<uint8_t> buffer;
        size_t position;
    };
    

    namespace implementation
    {
        template <typename T, template <typename> class R> T parse(R<BinaryBuffer>& reader);
        template <> GeometryCollection parse(Reader<BinaryBuffer>& reader);

        template <> BasePoint parse(Reader<BinaryBuffer>& reader) {
            double_t x = reader.read<double_t>();
            double_t y = reader.read<double_t>();
            return BasePoint(x, y);
        }

        template <> std::vector<BasePoint> parse(Reader<BinaryBuffer>& reader) {
            uint32_t pointCount = reader.read<uint32_t>();
            std::vector<BasePoint> ring;
            ring.reserve(pointCount);

            for (uint32_t i = 0; i < pointCount; ++i)
                ring.push_back(parse<BasePoint>(reader));

            return ring;
        }

        template <> Point parse(Reader<BinaryBuffer>& reader) {
            [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();

            assert(GeometryType::POINT == static_cast<GeometryType>(reader.read<uint32_t>()));

            return Point(parse<BasePoint>(reader));
        }

        template <> LineString parse(Reader<BinaryBuffer>& reader) {
            [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();

            assert(GeometryType::LINESTRING == static_cast<GeometryType>(reader.read<uint32_t>()));

            return LineString(parse<std::vector<BasePoint>>(reader));
        }

        template <> Polygon parse(Reader<BinaryBuffer>& reader) {
            [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();

            assert(GeometryType::POLYGON == static_cast<GeometryType>(reader.read<uint32_t>()));

            uint32_t ringCount = reader.read<uint32_t>();
            assert(ringCount > 0);

            std::vector<BasePoint> exteriorRing(parse<std::vector<BasePoint>>(reader));
            std::vector<std::vector<BasePoint>> interiorRings;

            for (uint32_t i = 1; i < ringCount; ++i) {
                interiorRings.push_back(parse<std::vector<BasePoint>>(reader));
            }

            return Polygon(std::move(exteriorRing), std::move(interiorRings));
        }

        template <> MultiPoint parse(Reader<BinaryBuffer>& reader) {
            [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();

            assert(GeometryType::MULTIPOINT == static_cast<GeometryType>(reader.read<uint32_t>()));

            uint32_t pointCount = reader.read<uint32_t>();
            std::vector<Point> points;
            points.reserve(pointCount);

            for (uint32_t i = 0; i < pointCount; ++i) {
                points.push_back(parse<Point>(reader));
            }

            return MultiPoint(std::move(points));
        }

        template <> MultiLineString parse(Reader<BinaryBuffer>& reader) {
            [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();

            assert(GeometryType::MULTILINESTRING == static_cast<GeometryType>(reader.read<uint32_t>()));

            uint32_t linestringCount = reader.read<uint32_t>();
            std::vector<LineString> linestrings;
            linestrings.reserve(linestringCount);

            for (uint32_t i = 0; i < linestringCount; ++i) {
                linestrings.push_back(parse<LineString>(reader));
            }

            return MultiLineString(std::move(linestrings));
        }

        template <> MultiPolygon parse(Reader<BinaryBuffer>& reader) {
            [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();

            assert(GeometryType::MULTIPOLYGON == static_cast<GeometryType>(reader.read<uint32_t>()));

            uint32_t polygonCount = reader.read<uint32_t>();
            std::vector<Polygon> polygons;
            polygons.reserve(polygonCount);

            for (uint32_t i = 0; i < polygonCount; ++i) {
                polygons.push_back(parse<Polygon>(reader));
            }

            return MultiPolygon(std::move(polygons));
        }

        template <> GenericGeometry parse(Reader<BinaryBuffer>& reader) {
            [[maybe_unused]] uint8_t endian = reader.peek<uint8_t>();
            GeometryType type  = static_cast<GeometryType>(reader.peek<uint32_t>(1));

            switch (type)
            {
            case GeometryType::POINT:
                return parse<Point>(reader);
            case GeometryType::LINESTRING:
                return parse<LineString>(reader);
            case GeometryType::POLYGON:
                return parse<Polygon>(reader);
            case GeometryType::MULTIPOINT:
                return parse<MultiPoint>(reader);
            case GeometryType::MULTILINESTRING:
                return parse<MultiLineString>(reader);
            case GeometryType::MULTIPOLYGON:
                return parse<MultiPolygon>(reader);
            case GeometryType::GEOMETRYCOLLECTION:
                return parse<GeometryCollection>(reader);
            default:
                throw std::runtime_error("Unkown geometry type " + std::to_string(static_cast<uint32_t>(type)));
            }
        }

        template <> GeometryCollection parse(Reader<BinaryBuffer>& reader) {
            [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();

            assert(GeometryType::GEOMETRYCOLLECTION == static_cast<GeometryType>(reader.read<uint32_t>()));

            uint32_t geometryCount = reader.read<uint32_t>();
            GeometryCollection geometries;
            geometries.reserve(geometryCount);

            for (uint32_t i = 0; i < geometryCount; ++i) {
                geometries.push_back(parse<GenericGeometry>(reader));
            }

            return geometries;
        }
    } // namespace implementation

    template <typename Geometry = GenericGeometry> Geometry from_wkb(BinaryBuffer& buffer) {
        Reader<BinaryBuffer> reader(buffer);

        return implementation::parse<Geometry>(reader);
    }

    BinaryBuffer to_wkb(const GenericGeometry& geometry) {

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

        BinaryBuffer buffer(std::visit(Y{compute_size}, geometry));
        size_t position = 0;

        auto write = [&buffer, &position]<typename T>(const T& value) -> void {
            assert(sizeof(T) + position <= buffer.size());

            memcpy(&(*buffer.begin()) + position, &value, sizeof(T));
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

        std::visit(Y{serialize}, geometry);

        return buffer;
    }

    
} // namespace tiledbsoma


#endif // TILEDBSOMA_GEOMETRY_IO_H