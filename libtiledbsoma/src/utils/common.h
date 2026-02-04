/**
 * @file   common.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This declares the common functions in the API
 */

#ifndef TILEDBSOMA_COMMON_H
#define TILEDBSOMA_COMMON_H

#include <functional>
#include <memory>
#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'
#include <string>
#include <string_view>
#include <tiledb/tiledb>

namespace tiledbsoma {

template <typename T>
using managed_unique_ptr = std::unique_ptr<T, std::function<void(T*)>>;

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

template <typename T>
class NoInitAlloc {
   public:
    using value_type = T;

    T* allocate(size_t n) {
        return static_cast<T*>(::operator new(n * sizeof(T)));
    }

    void deallocate(T* p, size_t n) noexcept {
        ::operator delete(p, sizeof(T) * n);
    }

    // Override construct to skip initialization
    template <typename U, typename... Args>
    void construct([[maybe_unused]] U* p, Args&&...) noexcept {
        // No initialization occurs
    }

    template <typename U>
    void destroy(U* p) noexcept {
        p->~U();
    }
};

template <typename T, typename U>
bool operator==(const NoInitAlloc<T>&, const NoInitAlloc<U>&) {
    return true;
}

template <typename T, typename U>
bool operator!=(const NoInitAlloc<T>&, const NoInitAlloc<U>&) {
    return false;
}

constexpr std::string_view SOMA_JOINID = "soma_joinid";

const std::string SOMA_OBJECT_TYPE_KEY = "soma_object_type";
const std::string ENCODING_VERSION_KEY = "soma_encoding_version";
const std::string ENCODING_VERSION_VAL = "1.1.0";
const std::string SPATIAL_ENCODING_VERSION_KEY = "soma_spatial_encoding_version";
const std::string SPATIAL_ENCODING_VERSION_VAL = "0.2.0";
const std::string SOMA_COORDINATE_SPACE_KEY = "soma_coordinate_space";
const std::string SOMA_GEOMETRY_COLUMN_NAME = "soma_geometry";
const std::string SOMA_GEOMETRY_DIMENSION_PREFIX = "tiledb__internal__";
const std::string ARROW_DATATYPE_METADATA_KEY = "dtype";

// SOMAColumn metadata keys
const std::string TILEDB_SOMA_SCHEMA_KEY = "tiledb_soma_schema";
const std::string TILEDB_SOMA_SCHEMA_VERSION = "0.0.1";
const std::string TILEDB_SOMA_SCHEMA_COL_KEY = "tiledb_columns";
const std::string TILEDB_SOMA_SCHEMA_COL_TYPE_KEY = "tiledb_column_type";
const std::string TILEDB_SOMA_SCHEMA_COL_DIM_KEY = "tiledb_dimensions";
const std::string TILEDB_SOMA_SCHEMA_COL_ATTR_KEY = "tiledb_attributes";

using MetadataValue = std::tuple<tiledb_datatype_t, uint32_t, const void*>;
enum MetadataInfo { dtype = 0, num, value };

using TimestampRange = std::pair<uint64_t, uint64_t>;

class TileDBSOMAError : public std::runtime_error {
   public:
    explicit TileDBSOMAError(const char* m)
        : std::runtime_error(m) {};
    explicit TileDBSOMAError(std::string m)
        : std::runtime_error(m.c_str()) {};

   public:
    virtual const char* what() const noexcept override {
        return std::runtime_error::what();
    }
};

// From
// https://github.com/TileDB-Inc/TileDB/blob/main/tiledb/common/scoped_executor.h
class ScopedExecutor final {
   public:
    /** Default constructor. */
    ScopedExecutor() = default;

    /**
     * Value constructor. Executes `fn` when this instance is
     * destructed.
     *
     * @param fn The function to execute on destruction.
     */
    explicit ScopedExecutor(std::function<void()>&& fn)
        : fn_(std::move(fn)) {
    }

    /** Move constructor. */
    ScopedExecutor(ScopedExecutor&& rhs) {
        fn_.swap(rhs.fn_);
    }

    /** Destructor. Executes `fn_`. */
    ~ScopedExecutor() {
        if (fn_) {
            fn_();
        }
    }

    ScopedExecutor(const ScopedExecutor&) = delete;

    ScopedExecutor& operator=(const ScopedExecutor&) = delete;

    // clang++ says:
    // error: constructor cannot be redeclared
    // ScopedExecutor(ScopedExecutor&&) = delete;

    ScopedExecutor& operator=(ScopedExecutor&&) = delete;

   private:
    /** The wrapped function to execute on destruction. */
    std::function<void()> fn_;
};

};  // namespace tiledbsoma

#endif  // TILEDBSOMA_COMMON_H
