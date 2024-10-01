#ifndef TILEDBSOMA_GEOMETRY_BASE_H
#define TILEDBSOMA_GEOMETRY_BASE_H

#include <vector>
#include <math.h>
#include <optional>
#include <memory>

namespace tiledbsoma {

struct BasePoint
{
    BasePoint(double_t x, double_t y, std::optional<double_t> z = std::nullopt, std::optional<double_t> m = std::nullopt) : x(x), y(y), z(z), m(m) {}

    double_t x;
    double_t y;
    std::optional<double_t> z;
    std::optional<double_t> m;
};
}

#endif  // TILEDBSOMA_GEOMETRY_BASE_H