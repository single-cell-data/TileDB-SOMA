#ifndef UTIL_H
#define UTIL_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'

#include <regex>
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

/**
 * @brief Check if the provided URI is a TileDB Cloud URI.
 *
 * @param uri URI to check
 * @return true URI is a TileBD Cloud URI
 * @return false URI is not a TileBD Cloud URI
 */
bool is_tiledb_uri(std::string_view uri);

/**
 * @brief Remove all trailing '/' from URI.
 *
 * @param uri URI
 * @return std::string URI without trailing '/'
 */
std::string rstrip_uri(std::string_view uri);

}  // namespace tiledbsc::util

#endif
