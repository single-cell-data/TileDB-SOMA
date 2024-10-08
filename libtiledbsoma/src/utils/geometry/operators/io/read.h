#ifndef TILEDBSOMA_GEOMETRY_READ_H
#define TILEDBSOMA_GEOMETRY_READ_H

#include <vector>
#include <string>
#include <cstring>
#include <variant>
#include <stdexcept>
#include <cassert>

#include "../../geometry.h"
#include "../../utils.h"

namespace tiledbsoma::geometry
{
    template <typename Input> 
    struct Reader {};
    
    template <> 
    struct Reader<BinaryBuffer>
    {
        Reader(const BinaryBuffer& buffer) : buffer(buffer), position(0) {}

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
        template <typename Geometry, template <typename> class R> Geometry parse(R<BinaryBuffer>& reader);

        template <> GeometryCollection parse(Reader<BinaryBuffer>& reader);
        template <> BasePoint parse(Reader<BinaryBuffer>& reader);
        template <> std::vector<BasePoint> parse(Reader<BinaryBuffer>& reader);
        template <> Point parse(Reader<BinaryBuffer>& reader);
        template <> LineString parse(Reader<BinaryBuffer>& reader);
        template <> Polygon parse(Reader<BinaryBuffer>& reader);
        template <> MultiPoint parse(Reader<BinaryBuffer>& reader);
        template <> MultiLineString parse(Reader<BinaryBuffer>& reader);
        template <> MultiPolygon parse(Reader<BinaryBuffer>& reader);
        template <> GenericGeometry parse(Reader<BinaryBuffer>& reader);
        template <> GeometryCollection parse(Reader<BinaryBuffer>& reader);
    } // namespace implementation

    template <typename Geometry = GenericGeometry> Geometry from_wkb(BinaryBuffer& buffer) {
        Reader<BinaryBuffer> reader(buffer);

        return implementation::parse<Geometry>(reader);
    }    
} // namespace tiledbsoma::geometry
#endif // TILEDBSOMA_GEOMETRY_READ_H