/**
 * @file   types.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 */

#ifndef COMMON_METADATA_TYPES_H
#define COMMON_METADATA_TYPES_H

#include <cstdint>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace tiledbsoma::common {
using MetadataSingleValueType =
    std::variant<int64_t, uint64_t, double, int32_t, uint32_t, float, int16_t, uint16_t, int8_t, uint8_t, bool>;

template <typename V, typename Seq>
struct MetadataValueTypeHelper;
template <typename V, int... S>
struct MetadataValueTypeHelper<V, std::integer_sequence<int, S...>> {
    using type = std::variant<
        std::vector<std::variant_alternative_t<S, V>>...,
        std::vector<std::byte>,
        std::string,
        std::variant_alternative_t<S, V>...>;
};

using MetadataValue = MetadataValueTypeHelper<MetadataSingleValueType, std::make_integer_sequence<int, std::variant_size_v<MetadataSingleValueType>>>::type;
}

#endif