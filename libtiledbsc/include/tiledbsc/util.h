#ifndef UTIL_H
#define UTIL_H

#include <tiledbsc/tiledbsc>

namespace tiledbsc::util {

using VarlenBufferPair =
    std::pair<std::vector<std::byte>, std::vector<uint64_t>>;

template <typename T>
VarlenBufferPair to_varlen_buffers(std::vector<T> data, bool arrow = true);

template <class T>
std::vector<T> to_vector(const std::span<T>& s) {
    return std::vector<T>(s.begin(), s.end());
}

}  // namespace tiledbsc::util

#endif
