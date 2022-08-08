#include "tiledbsc/soma.h"
#include "tiledbsc/logger_public.h"
#include "tiledbsc/util.h"

namespace tiledbsc {
using namespace tiledb;

//===================================================================
//= public static
//===================================================================

std::unique_ptr<SOMA> SOMA::open(
    std::string_view uri, std::shared_ptr<Context> ctx) {
    return std::make_unique<SOMA>(uri, ctx);
}

std::unique_ptr<SOMA> SOMA::open(std::string_view uri, const Config& config) {
    return std::make_unique<SOMA>(uri, std::make_shared<Context>(config));
}

//===================================================================
//= public non-static
//===================================================================

SOMA::SOMA(std::string_view uri, std::shared_ptr<Context> ctx)
    : ctx_(ctx)
    , uri_(util::rstrip_uri(uri)) {
}

std::unordered_map<std::string, std::string> SOMA::list_arrays() {
    // Allow only one thread to list the arrays
    std::lock_guard<std::mutex> lock(mtx_);

    if (array_uri_map_.empty()) {
        LOG_DEBUG(fmt::format("Listing arrays in SOMA '{}'", uri_));

        try {
            Group group(*ctx_, uri_, TILEDB_READ);
            build_uri_map(group);
        } catch (const std::exception& e) {
            throw TileDBSCError(fmt::format(
                "[SOMA] Error opening group URI='{}' : {}", uri_, e.what()));
        }
    }
    return array_uri_map_;
}

std::shared_ptr<Array> SOMA::open_array(const std::string& name) {
    // TODO: add option to open array without listing all arrays
    list_arrays();
    auto uri = array_uri_map_[name];
    LOG_DEBUG(fmt::format("Opening array '{}' from SOMA '{}'", name, uri_));

    try {
        return std::make_shared<Array>(*ctx_, uri, TILEDB_READ);
    } catch (const std::exception& e) {
        throw TileDBSCError(
            fmt::format("[SOMA] Error opening array '{}' : {}", uri, e.what()));
    }
}

//===================================================================
//= private non-static
//===================================================================

void SOMA::build_uri_map(Group& group, std::string_view parent) {
    // Iterate through all members in the group
    for (uint64_t i = 0; i < group.member_count(); i++) {
        auto member = group.member(i);
        auto path = parent.empty() ?
                        member.name().value() :
                        std::string(parent) + "/" + member.name().value();

        if (member.type() == Object::Type::Group) {
            // Member is a group, call recursively
            try {
                auto subgroup = Group(*ctx_, member.uri(), TILEDB_READ);
                build_uri_map(subgroup, path);
            } catch (const std::exception& e) {
                throw TileDBSCError(fmt::format(
                    "[SOMA] Error opening group URI='{}' : {}",
                    uri_,
                    e.what()));
            }
        } else {
            auto uri = member.uri();
            if (util::is_tiledb_uri(uri) && !util::is_tiledb_uri(uri_)) {
                // "Group member URI" is a TileDB Cloud URI, but the "SOMA
                // root URI" is *not* a TileDB Cloud URI. Build a "relative
                // group member URI"
                array_uri_map_[path] = uri_ + '/' + path;
                group_uri_override_ = true;
            } else {
                // Use the group member uri
                array_uri_map_[path] = uri;
            }
        }
    }
}
};  // namespace tiledbsc
