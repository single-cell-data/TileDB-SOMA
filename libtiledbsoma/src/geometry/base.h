#ifndef TILEDBSOMA_GEOMETRY_BASE_H
#define TILEDBSOMA_GEOMETRY_BASE_H

#include <math.h>
#include <memory>
#include <optional>
#include <vector>
#include <tiledbsoma_export.h>

namespace tiledbsoma::geometry {

struct TILEDBSOMA_EXPORT BasePoint {
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
