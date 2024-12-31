/**
 * @file   soma_dimension.h
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
 *   This file defines the SOMADimension class. SOMADimension acts as a wrapper
 * to a TileDB Dimension and implements function to perform queries as well as
 * core domain and current domain operations. It provides a common interface
 * identical to TileDB attributes and composite columns.
 */

#ifndef SOMA_DIMENSION_H
#define SOMA_DIMENSION_H

#include <algorithm>
#include <vector>

#include <tiledb/tiledb>
#include "soma_column.h"

namespace tiledbsoma {

using namespace tiledb;

class SOMADimension : public SOMAColumn {
   public:
    static std::shared_ptr<SOMADimension> create(
        std::shared_ptr<Context> ctx,
        ArrowSchema* schema,
        ArrowArray* array,
        const std::string& soma_type,
        std::string_view type_metadata,
        PlatformConfig platform_config);

    SOMADimension(Dimension dimension)
        : dimension(dimension) {
    }

    virtual inline std::string name() const override {
        return dimension.name();
    }

    virtual inline bool isIndexColumn() const override {
        return true;
    }

    virtual inline void select_columns(
        const std::unique_ptr<ManagedQuery>& query,
        bool if_not_empty = false) const override {
        query->select_columns(std::vector({dimension.name()}), if_not_empty);
    };

    virtual inline soma_column_datatype_t type() const override {
        return soma_column_datatype_t::SOMA_COLUMN_DIMENSION;
    }

    virtual inline std::optional<tiledb_datatype_t> domain_type()
        const override {
        return dimension.type();
    }

    virtual inline std::optional<tiledb_datatype_t> data_type() const override {
        return std::nullopt;
    }

    virtual inline std::optional<std::vector<Dimension>> tiledb_dimensions()
        override {
        return std::vector({dimension});
    }

    virtual inline std::optional<std::vector<Attribute>> tiledb_attributes()
        override {
        return std::nullopt;
    }

    virtual inline std::optional<std::vector<Enumeration>> tiledb_enumerations()
        override {
        return std::nullopt;
    }

    virtual ArrowArray* arrow_domain_slot(
        const SOMAContext& ctx,
        Array& array,
        enum Domainish kind) const override;

    virtual ArrowSchema* arrow_schema_slot(
        const SOMAContext& ctx, Array& array) override;

   protected:
    virtual void _set_dim_points(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::any& ranges) const override;

    virtual void _set_dim_ranges(
        const std::unique_ptr<ManagedQuery>& query,
        const SOMAContext& ctx,
        const std::any& ranges) const override;

    virtual void _set_current_domain_slot(
        NDRectangle& rectangle,
        std::span<const std::any> domain) const override;

    virtual std::pair<bool, std::string> _can_set_current_domain_slot(
        std::optional<NDRectangle>& rectangle,
        std::span<const std::any> new_domain) const override;

    virtual std::any _core_domain_slot() const override;

    virtual std::any _non_empty_domain_slot(Array& array) const override;

    virtual std::any _core_current_domain_slot(
        const SOMAContext& ctx, Array& array) const override;

    virtual std::any _core_current_domain_slot(
        NDRectangle& ndrect) const override;

   private:
    Dimension dimension;
};
}  // namespace tiledbsoma

#endif