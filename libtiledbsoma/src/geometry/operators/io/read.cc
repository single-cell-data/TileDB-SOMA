#include "read.h"

namespace tiledbsoma::geometry::implementation {
template <>
BasePoint parse(Reader<BinaryBuffer>& reader) {
    double_t x = reader.read<double_t>();
    double_t y = reader.read<double_t>();
    return BasePoint(x, y);
}

template <>
std::vector<BasePoint> parse(Reader<BinaryBuffer>& reader) {
    uint32_t pointCount = reader.read<uint32_t>();
    std::vector<BasePoint> ring;
    ring.reserve(pointCount);

    for (uint32_t i = 0; i < pointCount; ++i) {
        ring.push_back(parse<BasePoint>(reader));
    }

    return ring;
}

template <>
Point parse(Reader<BinaryBuffer>& reader) {
    [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();
    [[maybe_unused]] uint32_t type = reader.read<uint32_t>();
    assert(GeometryType::POINT == type);

    return Point(parse<BasePoint>(reader));
}

template <>
LineString parse(Reader<BinaryBuffer>& reader) {
    [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();
    [[maybe_unused]] uint32_t type = reader.read<uint32_t>();
    assert(GeometryType::LINESTRING == type);

    return LineString(parse<std::vector<BasePoint>>(reader));
}

template <>
Polygon parse(Reader<BinaryBuffer>& reader) {
    [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();
    [[maybe_unused]] uint32_t type = reader.read<uint32_t>();
    assert(GeometryType::POLYGON == type);

    uint32_t ringCount = reader.read<uint32_t>();
    assert(ringCount > 0);

    std::vector<BasePoint> exteriorRing(parse<std::vector<BasePoint>>(reader));
    std::vector<std::vector<BasePoint>> interiorRings;

    for (uint32_t i = 1; i < ringCount; ++i) {
        interiorRings.push_back(parse<std::vector<BasePoint>>(reader));
    }

    return Polygon(std::move(exteriorRing), std::move(interiorRings));
}

template <>
MultiPoint parse(Reader<BinaryBuffer>& reader) {
    [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();
    [[maybe_unused]] uint32_t type = reader.read<uint32_t>();
    assert(GeometryType::MULTIPOINT == type);

    uint32_t pointCount = reader.read<uint32_t>();
    std::vector<Point> points;
    points.reserve(pointCount);

    for (uint32_t i = 0; i < pointCount; ++i) {
        points.push_back(parse<Point>(reader));
    }

    return MultiPoint(std::move(points));
}

template <>
MultiLineString parse(Reader<BinaryBuffer>& reader) {
    [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();
    [[maybe_unused]] uint32_t type = reader.read<uint32_t>();
    assert(GeometryType::MULTILINESTRING == type);

    uint32_t linestringCount = reader.read<uint32_t>();
    std::vector<LineString> linestrings;
    linestrings.reserve(linestringCount);

    for (uint32_t i = 0; i < linestringCount; ++i) {
        linestrings.push_back(parse<LineString>(reader));
    }

    return MultiLineString(std::move(linestrings));
}

template <>
MultiPolygon parse(Reader<BinaryBuffer>& reader) {
    [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();
    [[maybe_unused]] uint32_t type = reader.read<uint32_t>();
    assert(GeometryType::MULTIPOLYGON == type);

    uint32_t polygonCount = reader.read<uint32_t>();
    std::vector<Polygon> polygons;
    polygons.reserve(polygonCount);

    for (uint32_t i = 0; i < polygonCount; ++i) {
        polygons.push_back(parse<Polygon>(reader));
    }

    return MultiPolygon(std::move(polygons));
}

template <>
GenericGeometry parse(Reader<BinaryBuffer>& reader) {
    [[maybe_unused]] uint8_t endian = reader.peek<uint8_t>();
    GeometryType type = static_cast<GeometryType>(reader.peek<uint32_t>(1));

    switch (type) {
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
            throw std::runtime_error(
                "Unkown geometry type " +
                std::to_string(static_cast<uint32_t>(type)));
    }
}

template <>
GeometryCollection parse(Reader<BinaryBuffer>& reader) {
    [[maybe_unused]] uint8_t endian = reader.read<uint8_t>();
    [[maybe_unused]] uint32_t type = reader.read<uint32_t>();
    assert(GeometryType::GEOMETRYCOLLECTION == type);

    uint32_t geometryCount = reader.read<uint32_t>();
    GeometryCollection geometries;
    geometries.reserve(geometryCount);

    for (uint32_t i = 0; i < geometryCount; ++i) {
        geometries.push_back(parse<GenericGeometry>(reader));
    }

    return geometries;
}
}  // namespace tiledbsoma::geometry::implementation