/**
 * @file   arrow.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *
 */

#ifndef ARROW_INTERFACE_H
#define ARROW_INTERFACE_H

#include <memory>
#include <string_view>
#include <tiledb/tiledb>

struct ArrowArray;
struct ArrowSchema;

namespace sparrow {
class array;
}

namespace tiledbsoma::arrow {

class Array {
   public:
    /* ********************************* */
    /*     CONSTRUCTORS & DESTRUCTORS    */
    /* ********************************* */

    /** No default constructor. Use the create factory method. */
    Array() = delete;

    /** Disable copy constructor. */
    Array(const Array&) = delete;
    /** Disable move constructor. */
    Array(Array&&) = delete;

    /** Destructor. */
    ~Array() = default;

    /* ********************************* */
    /*             OPERATORS             */
    /* ********************************* */

    /** Disable copy assignment. */
    Array& operator=(const Array&) = delete;
    /** Disable move assignment. */
    Array& operator=(const Array&&) = delete;

    /* ********************************* */
    /*                 API               */
    /* ********************************* */

    static std::shared_ptr<const Array> create(
        ArrowSchema* schema, ArrowArray* array);

    std::pair<ArrowArray, ArrowSchema> export_to_c() const noexcept;

    void export_to_c(ArrowArray* array, ArrowSchema* schema) const noexcept;

    std::string_view name() const noexcept;
    tiledb_datatype_t tdb_type() const noexcept;

   private:
    /**
     * Helper function to create shared_ptr<Enumeration> objects.
     * @tparam Args Argument types for the Enumeration constructor
     * @param args Arguments for the Enumeration constructor
     */
    template <typename... Args>
    static std::shared_ptr<Array> make_shared_array(Args&&... args) {
        struct shared_helper : public Array {
            shared_helper(Args&&... args)
                : Array(std::forward<Args>(args)...) {
            }
        };

        return std::make_shared<shared_helper>(std::forward<Args>(args)...);
    }

    Array(ArrowSchema* schema, ArrowArray* array);

    std::unique_ptr<sparrow::array> array_;
};
}  // namespace tiledbsoma::arrow

#endif