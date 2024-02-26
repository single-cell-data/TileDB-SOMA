#include <map>
#include <string>
#include <tiledb/tiledb>

#include "soma_array.h"
#include "soma_collection.h"
#include "soma_dataframe.h"
#include "soma_dense_ndarray.h"
#include "soma_experiment.h"
#include "soma_measurement.h"
#include "soma_sparse_ndarray.h"

namespace tiledbsoma {

using namespace tiledb;

std::unique_ptr<SOMAObject> SOMAObject::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp) {
    auto obj = tiledb::Object::object(*ctx->tiledb_ctx(), std::string(uri));

    if (obj.type() == tiledb::Object::Type::Array) {
        auto array_ = SOMAArray::open(
            mode, uri, ctx, "", {}, "auto", ResultOrder::automatic, timestamp);

        if (!array_->type().has_value())
            throw TileDBSOMAError("SOMAArray has no type info");

        if (array_->type() == "SOMADataFrame") {
            return std::make_unique<SOMADataFrame>(*array_);
        } else if (array_->type() == "SOMASparseNDArray") {
            return std::make_unique<SOMASparseNDArray>(*array_);
        } else if (array_->type() == "SOMADenseNDArray") {
            return std::make_unique<SOMADenseNDArray>(*array_);
        } else {
            throw TileDBSOMAError("Saw invalid SOMAArray type");
        }
    } else if (obj.type() == tiledb::Object::Type::Group) {
        auto group_ = SOMAGroup::open(mode, uri, ctx, "", timestamp);

        if (!group_->type().has_value())
            throw TileDBSOMAError("SOMAGroup has no type info");

        if (group_->type() == "SOMACollection") {
            return std::make_unique<SOMACollection>(*group_);
        } else if (group_->type() == "SOMAExperiment") {
            return std::make_unique<SOMAExperiment>(*group_);
        } else if (group_->type() == "SOMAMeasurement") {
            return std::make_unique<SOMAMeasurement>(*group_);
        } else {
            throw TileDBSOMAError("Saw invalid SOMAGroup type");
        }
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
