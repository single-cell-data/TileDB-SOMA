#include "multilinestring.h"

namespace tiledbsoma::geometry {
MultiLineString::MultiLineString(std::vector<LineString>&& linestrings)
    : linestrings(linestrings) {
}

MultiLineString::~MultiLineString() {
}
}  // namespace tiledbsoma::geometry
