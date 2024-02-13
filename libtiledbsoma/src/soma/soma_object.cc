#include <map>
#include <string>
#include <tiledb/tiledb>

#include "soma_array.h"
#include "soma_dataframe.h"
#include "soma_dense_ndarray.h"

namespace tiledbsoma {

using namespace tiledb;

std::unique_ptr<SOMAObject> SOMAObject::open(
    std::string_view uri,
    OpenMode mode,
    std::map<std::string, std::string> platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto ctx = std::make_shared<Context>(Config(platform_config));
    return SOMAObject::open(uri, mode, ctx, timestamp);
}

std::unique_ptr<SOMAObject> SOMAObject::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<Context> ctx,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto obj = tiledb::Object::object(*ctx, std::string(uri));

    if (obj.type() == tiledb::Object::Type::Array) {
        auto array_ = SOMAArray::open(
            mode, ctx, uri, "", {}, "auto", ResultOrder::automatic, timestamp);

        if (!array_->type().has_value())
            throw TileDBSOMAError(
                "Invalid SOMAObject passed to SOMAObject::open");

        if (*(array_->type()) == "SOMADataFrame")
            return std::make_unique<SOMADataFrame>(*array_);
        else if (*(array_->type()) == "SOMADenseNDArray")
            return std::make_unique<SOMADenseNDArray>(*array_);
        else
            throw TileDBSOMAError(
                "Invalid SOMAObject passed to SOMAObject::open");
    }

    throw TileDBSOMAError("Invalid TileDB object passed to SOMAObject::open");
}

const std::optional<std::string> SOMAObject::type() {
    auto soma_object_type = this->get_metadata("soma_object_type");

    if (!soma_object_type.has_value())
        return std::nullopt;

    const char* dtype = (const char*)std::get<MetadataInfo::value>(
        *soma_object_type);
    uint32_t sz = std::get<MetadataInfo::num>(*soma_object_type);

    return std::string(dtype, sz);
}

}  // namespace tiledbsoma
