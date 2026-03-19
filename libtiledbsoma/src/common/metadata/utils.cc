#include "utils.h"

#include <stdexcept>

#include "../logging/impl/logger.h"

namespace tiledbsoma::common {
MetadataValue decode_metadata(DataType type, uint32_t elements, const void* data) {
    auto decode = [&]<typename T>() -> MetadataValue {
        if constexpr (std::is_same_v<T, std::string>) {
            if (elements == 1 && data == nullptr) {
                return MetadataValue("");
            } else {
                return MetadataValue(std::string(static_cast<const char*>(data), elements));
            }
        } else if constexpr (std::is_same_v<T, bool>) {
            if (elements == 1) {
                return MetadataValue(reinterpret_cast<const uint8_t*>(data)[0] != 0);
            } else {
                return MetadataValue(
                    std::vector<bool>(
                        reinterpret_cast<const uint8_t*>(data), reinterpret_cast<const uint8_t*>(data) + elements));
            }
        } else {
            if (elements == 1) {
                return MetadataValue(reinterpret_cast<const T*>(data)[0]);
            } else {
                return MetadataValue(
                    std::vector<T>(reinterpret_cast<const T*>(data), reinterpret_cast<const T*>(data) + elements));
            }
        }
    };

    switch (type) {
        case common::DataType::boolean:
            return decode.template operator()<bool>();
        case common::DataType::int8:
            return decode.template operator()<int8_t>();
        case common::DataType::int16:
            return decode.template operator()<int16_t>();
        case common::DataType::int32:
            return decode.template operator()<int32_t>();
        case common::DataType::int64:
            return decode.template operator()<int64_t>();
        case common::DataType::uint8:
            return decode.template operator()<uint8_t>();
        case common::DataType::uint16:
            return decode.template operator()<uint16_t>();
        case common::DataType::uint32:
            return decode.template operator()<uint32_t>();
        case common::DataType::uint64:
            return decode.template operator()<uint64_t>();
        case common::DataType::float32:
            return decode.template operator()<float>();
        case common::DataType::float64:
            return decode.template operator()<double>();
        case common::DataType::string_ascii:
        case common::DataType::string_utf8:
            return decode.template operator()<std::string>();
        default:
            throw std::runtime_error(fmt::format("Unsupported metadata type '{}'", common::getName(type)));
    }
}
}  // namespace tiledbsoma::common