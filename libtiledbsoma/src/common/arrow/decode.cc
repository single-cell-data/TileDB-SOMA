#include "decode.h"

namespace tiledbsoma::common::arrow {
std::unordered_map<std::string, std::string> metadata_string_to_map(const char* metadata_str) {
    if (metadata_str == nullptr) {
        return std::unordered_map<std::string, std::string>();
    }

    // Current offset in the metadata pointer
    size_t offset = 0;
    std::unordered_map<std::string, std::string> metadata_map;

    int32_t element_count = reinterpret_cast<const int32_t*>(metadata_str)[0];
    offset += sizeof(int32_t);

    for (int32_t i = 0; i < element_count; ++i) {
        int32_t key_size = reinterpret_cast<const int32_t*>(metadata_str + offset)[0];
        offset += sizeof(int32_t);

        std::string key(metadata_str + offset, key_size);
        offset += key_size;

        int32_t value_size = reinterpret_cast<const int32_t*>(metadata_str + offset)[0];
        offset += sizeof(int32_t);

        std::string value(metadata_str + offset, value_size);
        offset += value_size;

        metadata_map.emplace(std::move(key), std::move(value));
    }

    return metadata_map;
}
}  // namespace tiledbsoma::common::arrow