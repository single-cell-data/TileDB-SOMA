/**
 * @file   soma_column.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2024 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAColumn class. SOMAColumn is an abstraction over
 * TileDB dimensions, attributes and combinations of them. It is designed to add
 * indexing capabilities to any datatype utilizing native TileDB dimensions
 * without exposing the internal indexing to the end user.
 */

#ifndef SOMA_COLUMN_H
#define SOMA_COLUMN_H

#include <any>
#include <format>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "enums.h"
#include "managed_query.h"
#include "nanoarrow/nanoarrow.hpp"
#include "utils/common.h"

namespace tiledbsoma {
using namespace tiledb;

class SOMAColumn {
   public:
    //===================================================================
    //= public non-static
    //===================================================================
    SOMAColumn() = default;
    SOMAColumn(const SOMAColumn&) = default;
    SOMAColumn(SOMAColumn&&) = default;
    SOMAColumn& operator=(const SOMAColumn&) = default;
    SOMAColumn& operator=(SOMAColumn&&) = default;

    virtual ~SOMAColumn() = default;

    /**
     * Get the SOMAColumn name as defined in schema.
     */
    virtual std::string name() const = 0;

    /**
     * If true, this column is used as index.
     *
     * @remark SOMAColumns used as indexes should define at least one TileDB
     * Dimension
     */
    virtual bool isIndexColumn() const = 0;

    /**
     * Get the TileDB Dimensions defined by the SOMAColumn object, if any.
     */
    virtual std::optional<std::vector<Dimension>> tiledb_dimensions() = 0;

    /**
     * Get the TileDB Attributes defined by the SOMAColumn object, if any.
     */
    virtual std::optional<std::vector<Attribute>> tiledb_attributes() = 0;

    /**
     * Get the TileDB Enumerations used by the SOMAColumn object, if any.
     */
    virtual std::optional<std::vector<Enumeration>> tiledb_enumerations() = 0;

    /**
     * Get the SOMAColumn type. Each subclass should define its own type.
     */
    virtual soma_column_datatype_t type() const = 0;

    /**
     * Get the datatype of the TileDB Dimensions if any. All dimensions must
     * have the same type.
     */
    virtual std::optional<tiledb_datatype_t> domain_type() const = 0;

    /**
     * Get the datatype of the TileDB Attributes if any. All attributes must
     * have the same type.
     */
    virtual std::optional<tiledb_datatype_t> data_type() const = 0;

    /**
     * @brief Select columns names to query (dim and attr). If the
     * `if_not_empty` parameter is `true`, the column will be selected iff the
     * list of selected columns is empty. This prevents a `select_columns` call
     * from changing an empty list (all columns) to a subset of columns.
     *
     * @param query the ManagedQuery object to modify
     * @param if_not_empty Prevent changing an "empty" selection of all columns
     */
    virtual void select_columns(
        const std::unique_ptr<ManagedQuery>& query,
        bool if_not_empty = false) const = 0;

    /**
     * Get the domain kind of the SOMAColumn as an ArrowArray for use with
     * R/Python API.
     */
    virtual ArrowArray* arrow_domain_slot(
        const SOMAContext& ctx,
        Array& array,
        enum Domainish which_kind) const = 0;

    /**
     * Get the SOMAColumn encoded as an ArrowSchema for use with R/Python API.
     */
    virtual ArrowSchema* arrow_schema_slot(
        const SOMAContext& ctx, Array& array) = 0;

    /**
     * Get the domain kind of the SOMAColumn.
     */
    template <typename T>
    std::pair<T, T> domain_slot(
        const SOMAContext& ctx, Array& array, enum Domainish which_kind) const {
        switch (which_kind) {
            case Domainish::kind_core_domain:
                return core_domain_slot<T>();
            case Domainish::kind_core_current_domain:
                return core_current_domain_slot<T>(ctx, array);
            case Domainish::kind_non_empty_domain:
                return non_empty_domain_slot<T>(array);
            default:
                throw std::runtime_error(
                    "internal coding error in SOMAArray::_core_domainish_slot: "
                    "unknown kind");
        }
    }

    /**
     * Set the current domain of this SOMAColumn.
     *
     * @param rectangle The current domain rectangle to modify.
     * @param domain A vector of void pointers to the the current domain data
     * buffers.
     */
    void set_current_domain_slot(
        NDRectangle& rectangle, const std::vector<const void*>& domain) const {
        if (!isIndexColumn()) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn] Column with name {} is not an index column",
                name()));
        }

        _set_current_domain_slot(rectangle, domain);
    }

    /**
     * @brief Set the dimension slice using one point
     *
     * @note Partitioning is not supported
     *
     * @tparam T
     * @param query
     * @param point
     */
    template <typename T>
    void set_dim_point(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const T& point) const {
        if (!isIndexColumn()) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn] Column with name {} is not an index column",
                name()));
        }

        T points[] = {point};
        this->_set_dim_points(
            query, ctx, std::make_any<std::span<T>>(std::span<T>(points)));
    }

    /**
     * @brief Set the dimension slice using multiple points
     *
     * @note Partitioning is not supported
     *
     * @tparam T
     * @param query
     * @param points
     */
    template <typename T>
    void set_dim_points(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        std::span<const T> points) const {
        if (!isIndexColumn()) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn] Column with name {} is not an index column",
                name()));
        }

        this->_set_dim_points(query, ctx, std::make_any<std::span<T>>(points));
    }

    /**
     * @brief Set the dimension slice using multiple ranges
     *
     * @note Partitioning is not supported
     *
     * @tparam T
     * @param query
     * @param ranges
     */
    template <typename T>
    void set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::vector<std::pair<T, T>>& ranges) const {
        if (!isIndexColumn()) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn] Column with name {} is not an index column",
                name()));
        }

        this->_set_dim_ranges(
            query, ctx, std::make_any<std::vector<std::pair<T, T>>>(ranges));
    }

    /**
     * Returns the core domain of this column.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     *   - soma maxdomain is core domain
     *   - core current domain is not accessed at the soma level
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    template <typename T>
    std::pair<T, T> core_domain_slot() const {
        try {
            return std::any_cast<std::pair<T, T>>(_core_domain_slot());
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        }
    }

    /**
     * Retrieves the non-empty domain from the array. This is the union of the
     * non-empty domains of the array fragments. Returns (0, 0) for empty
     * domains.
     */
    template <typename T>
    std::pair<T, T> non_empty_domain_slot(Array& array) const {
        try {
            return std::any_cast<std::pair<T, T>>(
                _non_empty_domain_slot(array));
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        };
    }

    /**
     * Returns the core current domain of this column.
     *
     * o For arrays with core current-domain support:
     *   - soma domain is core current domain
     *   - soma maxdomain is core domain
     * o For arrays without core current-domain support:
     *   - soma domain is core domain
     *   - soma maxdomain is core domain
     *   - core current domain is not accessed at the soma level
     *
     * @tparam T Domain datatype
     * @return Pair of [lower, upper] inclusive bounds.
     */
    template <typename T>
    std::pair<T, T> core_current_domain_slot(
        const SOMAContext& ctx, Array& array) const {
        try {
            return std::any_cast<std::pair<T, T>>(
                _core_current_domain_slot(ctx, array));
        } catch (const std::exception& e) {
            throw TileDBSOMAError(e.what());
        }
    }

   protected:
    virtual void _set_dim_points(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::any& points) const = 0;

    virtual void _set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::any& ranges) const = 0;

    virtual void _set_current_domain_slot(
        NDRectangle& rectangle,
        const std::vector<const void*>& domain) const = 0;

    virtual std::any _core_domain_slot() const = 0;

    virtual std::any _non_empty_domain_slot(Array& array) const = 0;

    virtual std::any _core_current_domain_slot(
        const SOMAContext& ctx, Array& array) const = 0;
};

template <>
std::pair<std::string, std::string> SOMAColumn::core_domain_slot<std::string>()
    const;

template <>
std::pair<std::string, std::string>
SOMAColumn::core_current_domain_slot<std::string>(
    const SOMAContext& ctx, Array& array) const;

}  // namespace tiledbsoma
#endif