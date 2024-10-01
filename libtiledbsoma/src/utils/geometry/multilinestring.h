#ifndef TILEDBSOMA_MULTILINESTRING_H
#define TILEDBSOMA_MULTILINESTRING_H

#include <vector>

#include "linestring.h"

namespace tiledbsoma
{
class MultiLineString {
public:
    MultiLineString(std::vector<LineString>&& linestrings);
    ~MultiLineString();

    std::vector<LineString> linestrings;
};
} // namespace tiledbsoma

#endif // TILEDBSOMA_MULTILINESTRING_H
