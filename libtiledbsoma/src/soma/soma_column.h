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
#include <vector>
#include <format>

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
     *
     * @param ctx
     * @param array
     * @param which_kind
     */
    virtual ArrowArray* arrow_domain_slot(
        const SOMAContext& ctx,
        Array& array,
        enum Domainish which_kind) const = 0;

    /**
     * Get the SOMAColumn encoded as an ArrowSchema for use with R/Python API.
     *
     * @param ctx
     * @param array
     */
    virtual ArrowSchema* arrow_schema_slot(
        const SOMAContext& ctx, Array& array) = 0;

    /**
     * Get the domain kind of the SOMAColumn.
     *
     * @tparam T
     * @param ctx
     * @param array
     * @param which_kind
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
     * @tparam T
     * @param rectangle The current domain rectangle to modify.
     * @param domain A vector of the n-dimensional domain in the form
     * [dim_0_min, dim_1_min, ..., dim_n_max]
     */
    template <typename T>
    void set_current_domain_slot(
        NDRectangle& rectangle, const std::vector<T>& domain) const {
        if (!isIndexColumn()) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_current_domain_slot] Column with name {} is "
                "not an index column",
                name()));
        }

        if (domain.size() % 2 != 0) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_current_domain_slot] Provided domain for "
                "column {} has missing values",
                name()));
        }

        std::vector<std::any> transformed_domain;
        size_t dim_count = domain.size() / 2;
        for (size_t i = 0; i < dim_count; ++i) {
            transformed_domain.push_back(std::make_any<std::array<T, 2>>(
                std::array<T, 2>({domain[i], domain[i + dim_count]})));
        }

        try {
            _set_current_domain_slot(rectangle, transformed_domain);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_current_domain_slot] Failed on \"{}\" with "
                "error \"{}\"",
                name(),
                e.what()));
        }
    }

    /**
     * Set the multi-type current domain of this SOMAColumn.
     *
     * @tparam T
     * @param rectangle The current domain rectangle to modify.
     * @param domain A vector holding std::arrays with 2 elements each [min,
     * max], casted as std::any
     */
    void set_current_domain_slot(
        NDRectangle& rectangle, const std::vector<std::any>& domain) const {
        if (!isIndexColumn()) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_current_domain_slot] Column with name {} is "
                "not an index column",
                name()));
        }

        try {
            _set_current_domain_slot(rectangle, domain);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_current_domain_slot] Failed on \"{}\" with "
                "error \"{}\"",
                name(),
                e.what()));
        }
    }

    /**
     * Test if the multi-type current domain of this SOMAColumn can be set with
     * the supplied new current domain.
     *
     * @tparam T
     * @param rectangle The current domain rectangle to modify.
     * @param domain A vector holding std::arrays with 2 elements each [min,
     * max], casted as std::any
     */
    std::pair<bool, std::string> can_set_current_domain_slot(
        std::optional<NDRectangle>& rectangle,
        const std::vector<std::any>& domain) const {
        if (!isIndexColumn()) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_current_domain_slot] Column with name {} is "
                "not an index column",
                name()));
        }

        try {
            return _can_set_current_domain_slot(rectangle, domain);
        } catch (const std::exception& e) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][can_set_current_domain_slot] Failed on \"{}\" "
                "with error \"{}\"",
                name(),
                e.what()));
        }
    }

    /**
     * @brief Set the dimension slice using one point
     *
     * @note Partitioning is not supported
     *
     * @tparam T
     * @param query
     * @param ctx
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

        try {
            this->_set_dim_points(
                query,
                ctx,
                std::make_any<std::span<const T>>(std::span<const T>(points)));
        } catch (const std::exception& e) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_dim_point] Failed on \"{}\" with error "
                "\"{}\"",
                name(),
                e.what()));
        }
    }

    /**
     * @brief Set the dimension slice using multiple points
     *
     * @note Partitioning is not supported
     *
     * @tparam T
     * @param query
     * @param ctx
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

        try {
            this->_set_dim_points(
                query, ctx, std::make_any<std::span<const T>>(points));
        } catch (const std::exception& e) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_dim_points] Failed on \"{}\" with error "
                "\"{}\"",
                name(),
                e.what()));
        }
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

        try {
            this->_set_dim_ranges(
                query,
                ctx,
                std::make_any<std::vector<std::pair<T, T>>>(ranges));
        } catch (const std::exception& e) {
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][set_dim_ranges] Failed on \"{}\" with error "
                "\"{}\"",
                name(),
                e.what()));
        }
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
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][core_domain_slot] Failed on \"{}\" with error "
                "\"{}\"",
                name(),
                e.what()));
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
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][non_empty_domain_slot] Failed on \"{}\" with "
                "error \"{}\"",
                name(),
                e.what()));
        }
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
            throw TileDBSOMAError(std::format(
                "[SOMAColumn][core_current_domain_slot] Failed on \"{}\" with "
                "error \"{}\"",
                name(),
                e.what()));
        }
    }

    /**
     * Returns the core current domain of this column from the supplied
     * NDRectangle.
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
    std::pair<T, T> core_current_domain_slot(NDRectangle& ndrect) const {
        try {
            return std::any_cast<std::pair<T, T>>(
                _core_current_domain_slot(ndrect));
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
        NDRectangle& rectangle, std::span<const std::any> domain) const = 0;

    virtual std::pair<bool, std::string> _can_set_current_domain_slot(
        std::optional<NDRectangle>& rectangle,
        std::span<const std::any> new_domain) const = 0;

    virtual std::any _core_domain_slot() const = 0;

    virtual std::any _non_empty_domain_slot(Array& array) const = 0;

    virtual std::any _core_current_domain_slot(
        const SOMAContext& ctx, Array& array) const = 0;

    virtual std::any _core_current_domain_slot(NDRectangle& ndrect) const = 0;
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