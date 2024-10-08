#ifndef TILEDBSOMA_GEOMETRY_WRITE_H
#define TILEDBSOMA_GEOMETRY_WRITE_H

#include <cassert>
#include <cstring>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "../../geometry.h"
#include "../../utils.h"

namespace tiledbsoma::geometry {
size_t wkb_size(const GenericGeometry& geometry);

void to_wkb(const GenericGeometry& geometry, uint8_t* buffer, size_t size);

BinaryBuffer to_wkb(const GenericGeometry& geometry);

}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_GEOMETRY_WRITE_H