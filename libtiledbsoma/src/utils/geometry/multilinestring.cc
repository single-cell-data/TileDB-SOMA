#include "multilinestring.h"

namespace tiledbsoma
{
    MultiLineString::MultiLineString(std::vector<LineString>&& linestrings) : linestrings(linestrings) {}

    MultiLineString::~MultiLineString() {}
} // namespace tiledbsoma
