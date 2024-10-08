#ifndef TILEDBSOMA_GEOMETRY_UTILS_H
#define TILEDBSOMA_GEOMETRY_UTILS_H

#include <utility>

namespace tiledbsoma::geometry {

template <typename... Base>
struct Operator : Base... {
    using Base::operator()...;
};

template <typename... T>
Operator(T...) -> Operator<T...>;

template <typename G>
struct Recurse {
    template <typename... X>
    decltype(auto) operator()(X&&... x) const& {
        return g(*this, std::forward<X>(x)...);
    }
    G g;
};

template <typename G>
Recurse(G) -> Recurse<G>;

}  // namespace tiledbsoma::geometry

#endif