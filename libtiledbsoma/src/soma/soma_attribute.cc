/**
 * @file   soma_attribute.cc
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
 *   This file defines the SOMAAttribute class.
 */

#include "soma_attribute.h"

namespace tiledbsoma {
std::shared_ptr<SOMAAttribute> SOMAAttribute::create(
    std::shared_ptr<Context> ctx,
    ArrowSchema* schema,
    std::string_view type_metadata,
    PlatformConfig platform_config) {
    auto attribute = ArrowAdapter::tiledb_attribute_from_arrow_schema(
        ctx, schema, type_metadata, platform_config);

    return std::make_shared<SOMAAttribute>(
        SOMAAttribute(attribute.first, attribute.second));
}

void SOMAAttribute::_set_dim_points(
    const std::unique_ptr<ManagedQuery>&,
    const SOMAContext&,
    const std::any&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_set_dim_points] Column with name {} is not an index "
        "column",
        name()));
}

void SOMAAttribute::_set_dim_ranges(
    const std::unique_ptr<ManagedQuery>&,
    const SOMAContext&,
    const std::any&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_set_dim_ranges] Column with name {} is not an index "
        "column",
        name()));
}

void SOMAAttribute::_set_current_domain_slot(
    NDRectangle&, std::span<const std::any>) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_set_current_domain_slot] Column with name {} is not "
        "an index column",
        name()));
}

std::pair<bool, std::string> SOMAAttribute::_can_set_current_domain_slot(
    std::optional<NDRectangle>&, std::span<const std::any>) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_set_current_domain_slot] Column with name {} is not "
        "an index column",
        name()));
};

std::any SOMAAttribute::_core_domain_slot() const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_core_domain_slot] Column with name {} is not an "
        "index column",
        name()));
}

std::any SOMAAttribute::_non_empty_domain_slot(Array&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_non_empty_domain_slot] Column with name {} is not an "
        "index column",
        name()));
}

std::any SOMAAttribute::_core_current_domain_slot(
    const SOMAContext&, Array&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_core_current_domain_slot] Column with name {} is not "
        "an index column",
        name()));
}

std::any SOMAAttribute::_core_current_domain_slot(NDRectangle&) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][_core_current_domain_slot] Column with name {} is not "
        "an index column",
        name()));
}

ArrowArray* SOMAAttribute::arrow_domain_slot(
    const SOMAContext&, Array&, enum Domainish) const {
    throw TileDBSOMAError(std::format(
        "[SOMAAttribute][arrow_domain_slot] Column with name {} is not an "
        "index column",
        name()));
}

ArrowSchema* SOMAAttribute::arrow_schema_slot(
    const SOMAContext& ctx, Array& array) {
    return ArrowAdapter::arrow_schema_from_tiledb_attribute(
               attribute, *ctx.tiledb_ctx(), array)
        .release();
}
}  // namespace tiledbsoma