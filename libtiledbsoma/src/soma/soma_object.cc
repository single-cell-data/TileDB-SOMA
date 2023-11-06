#include <map>
#include <string>
#include <tiledb/tiledb>

#include "soma_array.h"
#include "soma_dataframe.h"

namespace tiledbsoma {

using namespace tiledb;

std::unique_ptr<SOMAObject> SOMAObject::open(
    std::string uri,
    OpenMode mode,
    std::map<std::string, std::string> platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto ctx = std::make_shared<Context>(Config(platform_config));
    return SOMAObject::open(uri, mode, ctx, timestamp);
}

std::unique_ptr<SOMAObject> SOMAObject::open(
    std::string uri,
    OpenMode mode,
    std::shared_ptr<Context> ctx,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto obj = tiledb::Object::object(*ctx, uri);

    if (obj.type() == tiledb::Object::Type::Array) {
        auto array_ = SOMAArray::open(
            mode, ctx, uri, "", {}, "auto", ResultOrder::automatic, timestamp);

        if (array_->type() == "SOMADataFrame")
            return std::make_unique<SOMADataFrame>(std::move(array_));
        else
            throw TileDBSOMAError(
                "Invalid SOMAObject passed to SOMAObject::open");
    }

    throw TileDBSOMAError("Invalid TileDB object passed to SOMAObject::open");
}

}  // namespace tiledbsoma