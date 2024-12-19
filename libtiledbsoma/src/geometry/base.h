#ifndef TILEDBSOMA_GEOMETRY_BASE_H
#define TILEDBSOMA_GEOMETRY_BASE_H

#include <math.h>
#include <memory>
#include <optional>
#include <vector>

namespace tiledbsoma::geometry {

struct BasePoint {
    BasePoint(
        double_t x,
        double_t y,
        std::optional<double_t> z = std::nullopt,
        std::optional<double_t> m = std::nullopt)
        : x(x)
        , y(y)
        , z(z)
        , m(m) {
    }

    virtual ~BasePoint() = default;

    double_t x;
    double_t y;
    std::optional<double_t> z;
    std::optional<double_t> m;
};
}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_GEOMETRY_BASE_H
