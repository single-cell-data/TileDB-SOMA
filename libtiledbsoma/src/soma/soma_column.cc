/**
 * @file   soma_column.cc
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAColumn class.
 */

#include "soma_column.h"

namespace tiledbsoma {

template <>
std::pair<std::string, std::string> SOMAColumn::core_domain_slot<std::string>()
    const {
    return std::pair<std::string, std::string>("", "");
}

template <>
std::pair<std::string, std::string>
SOMAColumn::core_current_domain_slot<std::string>(
    const SOMAContext& ctx, Array& array) const {
    // Here is an intersection of a few oddities:
    //
    // * Core domain for string dims must be a nullptr pair; it cannot
    // be
    //   anything else.
    // * TileDB-Py shows this by using an empty-string pair, which we
    //   imitate.
    // * Core current domain for string dims must _not_ be a nullptr
    // pair.
    // * In TileDB-SOMA, unless the user specifies otherwise, we use ""
    // for
    //   min and "\x7f" for max. (We could use "\x7f" but that causes
    //   display problems in Python.)
    //
    // To work with all these factors, if the current domain is the
    // default
    // "" to "\x7f", return an empty-string pair just as we do for
    // domain. (There was some pre-1.15 software using "\xff" and it's
    // super-cheap to check for that as well.)
    try {
        std::pair<std::string, std::string>
            current_domain = std::any_cast<std::pair<std::string, std::string>>(
                _core_current_domain_slot(ctx, array));

        if (current_domain.first == "" && (current_domain.second == "\x7f" ||
                                           current_domain.second == "\xff")) {
            return std::pair<std::string, std::string>("", "");
        }

        return current_domain;
    } catch (const std::exception& e) {
        throw TileDBSOMAError(e.what());
    }
}
}  // namespace tiledbsoma
