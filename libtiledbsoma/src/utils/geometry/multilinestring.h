#ifndef TILEDBSOMA_MULTILINESTRING_H
#define TILEDBSOMA_MULTILINESTRING_H

#include <vector>

#include "linestring.h"

namespace tiledbsoma::geometry {
class MultiLineString {
   public:
    MultiLineString(
        std::vector<LineString>&& linestring = std::vector<LineString>());
    ~MultiLineString();

    std::vector<LineString> linestrings;
};
}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_MULTILINESTRING_H
