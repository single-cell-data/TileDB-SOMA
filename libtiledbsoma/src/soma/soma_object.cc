#include <map>
#include <string>
#include <tiledb/tiledb>

#include "../utils/util.h"
#include "common/logging/impl/logger.h"
#include "soma_array.h"
#include "soma_collection.h"
#include "soma_dataframe.h"
#include "soma_dense_ndarray.h"
#include "soma_experiment.h"
#include "soma_geometry_dataframe.h"
#include "soma_measurement.h"
#include "soma_multiscale_image.h"
#include "soma_point_cloud_dataframe.h"
#include "soma_scene.h"
#include "soma_sparse_ndarray.h"

namespace tiledbsoma {

using namespace tiledb;

std::unique_ptr<SOMAObject> SOMAObject::open(
    std::string_view uri,
    OpenMode mode,
    std::shared_ptr<SOMAContext> ctx,
    std::optional<TimestampRange> timestamp,
    std::optional<std::string> soma_type) {
    auto tiledb_type = Object::Type::Invalid;
    if (soma_type.has_value()) {
        if (soma_type == "SOMAArray") {
            tiledb_type = Object::Type::Array;
        } else if (soma_type == "SOMAGroup") {
            tiledb_type = Object::Type::Group;
        } else {
            throw TileDBSOMAError(fmt::format("Internal error: Invalid soma base type '{}'.", soma_type.value()));
        }
    } else {
        tiledb_type = Object::object(*ctx->tiledb_ctx(), std::string(uri)).type();
    }

    switch (tiledb_type) {
        case Object::Type::Array: {
            auto array_ = SOMAArray::open(mode, uri, ctx, timestamp);
            auto array_type = array_->type();
            if (!array_type.has_value()) {
                throw TileDBSOMAError(
                    fmt::format(
                        "Unable to open the TileDB array at '{}'. The array is missing the required '{}' metadata key.",
                        uri,
                        SOMA_OBJECT_TYPE_KEY));
            }
            if (array_type == "SOMADataFrame") {
                array_->check_encoding_version();
                return std::make_unique<SOMADataFrame>(*array_);
            } else if (array_type == "SOMASparseNDArray") {
                array_->check_encoding_version();
                return std::make_unique<SOMASparseNDArray>(*array_);
            } else if (array_type == "SOMADenseNDArray") {
                array_->check_encoding_version();
                return std::make_unique<SOMADenseNDArray>(*array_);
            } else if (array_type == "SOMAPointCloudDataFrame") {
                array_->check_encoding_version();
                return std::make_unique<SOMAPointCloudDataFrame>(*array_);
            } else if (array_type == "SOMAGeometryDataFrame") {
                array_->check_encoding_version();
                return std::make_unique<SOMAGeometryDataFrame>(*array_);
            } else {
                throw TileDBSOMAError(
                    fmt::format(
                        "Unable to open the TileDB array at '{}' with unrecognized SOMA array type '{}'.",
                        uri,
                        array_->type().value()));
            }
        }
        case Object::Type::Group: {
            auto group_ = SOMAGroup::open(mode, uri, ctx, "", timestamp);
            auto group_type = group_->type();
            if (!group_type.has_value()) {
                throw TileDBSOMAError(
                    fmt::format(
                        "Unable to open the TileDB group at '{}'. The group is missing the required '{}' metadata "
                        "key.",
                        uri,
                        SOMA_OBJECT_TYPE_KEY));
            }
            if (group_type == "SOMACollection") {
                group_->check_encoding_version();
                return std::make_unique<SOMACollection>(*group_);
            } else if (group_type == "SOMAExperiment") {
                group_->check_encoding_version();
                return std::make_unique<SOMAExperiment>(*group_);
            } else if (group_type == "SOMAMeasurement") {
                group_->check_encoding_version();
                return std::make_unique<SOMAMeasurement>(*group_);
            } else if (group_type == "SOMAScene") {
                group_->check_encoding_version();
                return std::make_unique<SOMAScene>(*group_);
            } else if (group_type == "SOMAMultiscaleImage") {
                group_->check_encoding_version();
                return std::make_unique<SOMAMultiscaleImage>(*group_);
            } else {
                throw TileDBSOMAError(
                    fmt::format(
                        "Unable to open the TileDB group at '{}' with unrecognized SOMA group type '{}'.",
                        uri,
                        group_->type().value()));
            }
        }
        default:
            throw TileDBSOMAError(fmt::format("The object at URI '{}' is not a valid TileDB array or group.", uri));
    }
}

const std::optional<std::string> SOMAObject::type() {
    auto soma_object_type = get_metadata(SOMA_OBJECT_TYPE_KEY);

    if (!soma_object_type.has_value()) {
        return std::nullopt;
    }

    const char* dtype = static_cast<const char*>(std::get<MetadataInfo::value>(*soma_object_type));
    uint32_t sz = std::get<MetadataInfo::num>(*soma_object_type);

    return std::string(dtype, sz);
}

const std::optional<std::string> SOMAObject::encoding_version() {
    auto encoding_version = get_metadata(ENCODING_VERSION_KEY);

    if (!encoding_version.has_value()) {
        return std::nullopt;
    }

    const char* dtype = static_cast<const char*>(std::get<MetadataInfo::value>(*encoding_version));
    uint32_t sz = std::get<MetadataInfo::num>(*encoding_version);

    return std::string(dtype, sz);
}

void SOMAObject::check_encoding_version() {
    auto version = encoding_version();
    if (!version.has_value()) {
        throw TileDBSOMAError(
            fmt::format(
                "Unable to open object at '{}'. Object is missing required '{}' metadata key.",
                std::string(uri()),
                ENCODING_VERSION_KEY));
    }
    if (version.value() != "1" && version.value() != "1.1.0") {
        throw TileDBSOMAError(
            fmt::format(
                "Unable to open object at '{}' with unsupported encoding version '{}'. Try updating "
                "to a more recent release of TileDB-SOMA.",
                std::string(uri()),
                version.value()));
    }
};

tiledb_object_t SOMAObject::tiledb_type_from_soma_type(const std::string& soma_type) {
    const std::map<std::string, tiledb_object_t> typeMap = {
        {"SOMAArray", tiledb_object_t::TILEDB_ARRAY},
        {"SOMACollection", tiledb_object_t::TILEDB_GROUP},
        {"SOMADataFrame", tiledb_object_t::TILEDB_ARRAY},
        {"SOMADenseNDArray", tiledb_object_t::TILEDB_ARRAY},
        {"SOMAExperiment", tiledb_object_t::TILEDB_GROUP},
        {"SOMAGeometryDataFrame", tiledb_object_t::TILEDB_ARRAY},
        {"SOMAGroup", tiledb_object_t::TILEDB_GROUP},
        {"SOMAMeasurement", tiledb_object_t::TILEDB_GROUP},
        {"SOMAMultiscaleImage", tiledb_object_t::TILEDB_GROUP},
        {"SOMAPointCloudDataFrame", tiledb_object_t::TILEDB_ARRAY},
        {"SOMAScene", tiledb_object_t::TILEDB_GROUP},
        {"SOMASparseNDArray", tiledb_object_t::TILEDB_ARRAY},
    };
    const std::map<std::string, tiledb_object_t>::const_iterator iTileDBType = typeMap.find(soma_type);
    if (iTileDBType == typeMap.end())
        return tiledb_object_t::TILEDB_INVALID;
    return iTileDBType->second;
};

}  // namespace tiledbsoma
