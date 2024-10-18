/**
 * @file   util.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file defines utilities.
 */

#include "utils/util.h"
#include <cstring>
#include "../geometry/geometry.h"
#include "../geometry/operators/envelope.h"
#include "../geometry/operators/io/write.h"
#include "logger.h"

namespace tiledbsoma::util {

template <typename T>
VarlenBufferPair to_varlen_buffers(std::vector<T> data, bool arrow) {
    size_t nbytes = 0;
    for (auto& elem : data) {
        nbytes += elem.size();
    }

    std::string result;
    std::vector<uint64_t> offsets(data.size() + 1);
    size_t offset = 0;
    size_t idx = 0;

    for (auto& elem : data) {
        result += elem;
        offsets[idx++] = offset;
        offset += elem.size();
    }
    offsets[idx] = offset;

    // Remove extra arrow offset when creating buffers for TileDB write
    if (!arrow) {
        offsets.pop_back();
    }

    return {result, offsets};
}

template VarlenBufferPair to_varlen_buffers(
    std::vector<std::string>, bool arrow);

bool is_tiledb_uri(std::string_view uri) {
    return uri.find("tiledb://") == 0;
}

std::string rstrip_uri(std::string_view uri) {
    return std::regex_replace(std::string(uri), std::regex("/+$"), "");
}

std::vector<uint8_t> cast_bit_to_uint8(ArrowSchema* schema, ArrowArray* array) {
    if (strcmp(schema->format, "b") != 0) {
        throw TileDBSOMAError(fmt::format(
            "_cast_bit_to_uint8 expected column format to be 'b' but saw {}",
            schema->format));
    }

    const void* data;
    if (array->n_buffers == 3) {
        data = array->buffers[2];
    } else {
        data = array->buffers[1];
    }

    std::vector<uint8_t> casted;
    for (int64_t i = 0; i * 8 < array->length; ++i) {
        uint8_t byte = ((uint8_t*)data)[i];
        for (int64_t j = 0; j < 8; ++j) {
            casted.push_back((uint8_t)((byte >> j) & 0x01));
        }
    }
    return casted;
}

std::vector<ArrowArray*> cast_vertices_to_wkb(
    ArrowArray* array, std::vector<std::string> spatial_axes) {
    // Initialize a vector to hold all the Arrow tables containing the
    // transformed geometry data
    ArrowError error;
    std::vector<ArrowArray*> arrays({(ArrowArray*)malloc(sizeof(ArrowArray))});
    NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
        arrays[0], ArrowType::NANOARROW_TYPE_LARGE_BINARY));

    for (auto axis : spatial_axes) {
        // Min spatial axis
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            arrays.emplace_back((ArrowArray*)malloc(sizeof(ArrowArray))),
            ArrowType::NANOARROW_TYPE_DOUBLE));

        // Max spatial axis
        NANOARROW_THROW_NOT_OK(ArrowArrayInitFromType(
            arrays.emplace_back((ArrowArray*)malloc(sizeof(ArrowArray))),
            ArrowType::NANOARROW_TYPE_DOUBLE));
    }

    // Large list of doubles
    const uint32_t* offset = static_cast<const uint32_t*>(array->buffers[1]);
    const double_t* data = static_cast<const double_t*>(
        array->children[0]->buffers[1]);

    size_t wkb_buffer_size = 0;
    std::vector<geometry::GenericGeometry> geometries;

    for (int64_t index = 0; index < array->length; ++index) {
        int64_t stop_index = index < array->length - 1 ?
                                 offset[index + 1] :
                                 array->children[0]->length;

        std::vector<geometry::BasePoint> ring;
        for (int64_t j = offset[index]; j < stop_index; j += 2) {
            ring.push_back(geometry::BasePoint(data[j], data[j + 1]));
        }

        geometries.push_back(
            geometry::GenericGeometry(geometry::Polygon(std::move(ring))));
        wkb_buffer_size += wkb_size(geometries.back());
    }

    NANOARROW_THROW_NOT_OK(ArrowArrayReserve(arrays.front(), wkb_buffer_size));
    NANOARROW_THROW_NOT_OK(ArrowArrayStartAppending(arrays.front()));
    for (size_t i = 1; i < arrays.size(); ++i) {
        NANOARROW_THROW_NOT_OK(ArrowArrayReserve(arrays[i], array->length));
        NANOARROW_THROW_NOT_OK(ArrowArrayStartAppending(arrays[i]));
    }

    for (auto& geometry : geometries) {
        geometry::BinaryBuffer wkb = geometry::to_wkb(geometry);
        geometry::Envelope envelope = geometry::envelope(geometry);

        ArrowBufferView wkb_view;
        wkb_view.data.data = wkb.data();
        wkb_view.size_bytes = (int64_t)wkb.size();

        NANOARROW_THROW_NOT_OK(ArrowArrayAppendBytes(arrays.front(), wkb_view));

        for (size_t i = 0; i < spatial_axes.size(); ++i) {
            NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(
                arrays[2 * i + 1], envelope.range.at(i).first));
            NANOARROW_THROW_NOT_OK(ArrowArrayAppendDouble(
                arrays[2 * i + 2], envelope.range.at(i).second));
        }
    }

    NANOARROW_THROW_NOT_OK(
        ArrowArrayFinishBuildingDefault(arrays.front(), &error));
    for (size_t i = 1; i < arrays.size(); ++i) {
        NANOARROW_THROW_NOT_OK(
            ArrowArrayFinishBuildingDefault(arrays[i], &error));
    }

    return arrays;
}

};  // namespace tiledbsoma::util
