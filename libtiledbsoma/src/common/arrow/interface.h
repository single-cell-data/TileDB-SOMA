/**
 * @file   interface.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 */

#ifndef COMMON_ARROW_INTERFACE_H
#define COMMON_ARROW_INTERFACE_H

#include <span>

namespace tiledbsoma::common::arrow {
class IArrowBufferStorage {
   public:
    virtual ~IArrowBufferStorage() {
    }

    virtual std::span<const std::byte> data() const = 0;
    virtual std::span<const std::byte> offsets() const = 0;
    virtual std::span<const std::byte> validity() const = 0;

    size_t length() const;
    size_t null_count() const;

   protected:
    size_t length_;
    size_t null_count_;
    size_t data_size_;
    size_t offsets_size_;
};
}  // namespace tiledbsoma::common::arrow

#endif