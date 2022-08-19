/**
 * @file   soma_collection.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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
 * This file defines the SOMA collection class.
 */

#include <regex>

#include "tiledbsc/soma_collection.h"
#include "tiledbsc/util.h"

namespace tiledbsc {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMACollection> SOMACollection::open(
    std::string_view uri, std::shared_ptr<Context> ctx) {
    return std::make_unique<SOMACollection>(uri, ctx);
}

std::unique_ptr<SOMACollection> SOMACollection::open(
    std::string_view uri, const Config& config) {
    return std::make_unique<SOMACollection>(
        uri, std::make_shared<Context>(config));
}

//===================================================================
//= public non-static
//===================================================================

SOMACollection::SOMACollection(
    std::string_view uri, std::shared_ptr<Context> ctx)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(uri)) {
}

std::unordered_map<std::string, std::string> SOMACollection::list_somas() {
    if (soma_uri_map_.empty()) {
        Group group(*ctx_, uri_, TILEDB_READ);
        build_uri_map(group);
    }
    return soma_uri_map_;
}

std::unordered_map<std::string, std::shared_ptr<SOMA>>
SOMACollection::get_somas() {
    if (soma_map_.empty()) {
        for (auto& [name, uri] : list_somas()) {
            soma_map_[name] = SOMA::open(uri, ctx_);
        }
    }

    return soma_map_;
}

//===================================================================
//= private non-static
//===================================================================

void SOMACollection::build_uri_map(Group& group, std::string_view parent) {
    // Iterate through all members in the group
    for (uint64_t i = 0; i < group.member_count(); i++) {
        auto member = group.member(i);
        auto path = parent.empty() ?
                        member.name().value() :
                        std::string(parent) + "/" + member.name().value();

        if (member.type() == Object::Type::Group) {
            auto subgroup = Group(*ctx_, member.uri(), TILEDB_READ);
            // Determine if the subgroup is a SOMA or a nested SOCO

            // Read group metadata "__soma_object_type__"
            tiledb_datatype_t value_type;
            uint32_t value_num;
            const void* value;
            subgroup.get_metadata(
                "__soma_object_type__", &value_type, &value_num, &value);

            bool subgroup_is_soma;
            if (value) {
                // Found metadata, check if object_type is "SOMA"
                std::string_view object_type(
                    static_cast<const char*>(value), value_num);
                subgroup_is_soma = object_type == "SOMA";
            } else {
                // Metadata not found, infer the group is a SOMA if it contains
                // an obs array
                subgroup_is_soma = subgroup.member("obs").type() ==
                                   Object::Type::Array;
            }

            if (subgroup_is_soma) {
                // Since SOCO members may reference tiledb:// or file://
                // SOMAs outside of the SOCO relative directory structure,
                // do not modify the group member URI like the SOMA class does.
                soma_uri_map_[path] = member.uri();
            } else {
                // Member is a SOCO, call recursively
                build_uri_map(subgroup);
            }
        }
    }
}
};  // namespace tiledbsc
