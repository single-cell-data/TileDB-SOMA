#include <regex>

#include "tiledbsc/soma_collection.h"

namespace tiledbsc {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

SOMACollection SOMACollection::open(std::string_view uri, Context ctx) {
    return SOMACollection(uri, ctx);
}

//===================================================================
//= public non-static
//===================================================================

SOMACollection::SOMACollection(std::string_view uri, Context ctx)
    : ctx_(ctx) {
    // Remove all trailing /
    // TODO: move this to utils
    uri_ = std::regex_replace(std::string(uri), std::regex("/+$"), "");
}

std::unordered_map<std::string, std::string> SOMACollection::list_somas() {
    if (soma_uri_map_.empty()) {
        Group group(ctx_, uri_, TILEDB_READ);
        build_uri_map(group);
    }
    return soma_uri_map_;
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
            auto subgroup = Group(ctx_, member.uri(), TILEDB_READ);
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
