#ifndef COMMON_CONCEPTS_H
#define COMMON_CONCEPTS_H

#include <concepts>
#include <memory>

template <typename T>
concept is_data_buffer = std::same_as<std::unique_ptr<std::byte[]>, T> ||
                         (std::is_pointer_v<T> &&
                          (std::same_as<std::remove_const_t<std::remove_pointer_t<T>>, void> ||
                           std::same_as<std::remove_const_t<std::remove_pointer_t<T>>, std::byte> ||
                           std::integral<std::remove_const_t<std::remove_pointer_t<T>>> ||
                           std::floating_point<std::remove_const_t<std::remove_pointer_t<T>>>));

template <typename T>
concept is_offset_buffer = std::same_as<T, std::unique_ptr<uint64_t[]>> ||
                           (std::is_pointer_v<T> &&
                            (std::same_as<std::remove_const_t<std::remove_pointer_t<T>>, uint32_t> ||
                             std::same_as<std::remove_const_t<std::remove_pointer_t<T>>, uint64_t>));

#endif