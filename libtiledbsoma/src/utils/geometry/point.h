#ifndef TILEDBSOMA_POINT_H
#define TILEDBSOMA_POINT_H

#include <math.h>
#include <optional>

#include "base.h"

namespace tiledbsoma::geometry {

class Point : public BasePoint {
   public:
    Point();
    Point(
        double_t x,
        double_t y,
        std::optional<double_t> z = std::nullopt,
        std::optional<double_t> m = std::nullopt);
    Point(BasePoint&& point);

    ~Point();
};
}  // namespace tiledbsoma::geometry

#endif  // TILEDBSOMA_POINT_H
