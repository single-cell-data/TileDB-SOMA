#include <map>
#include <string>
#include <tiledb/tiledb>

#include "../utils/logger.h"
#include "../utils/util.h"
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
    if (soma_type == std::nullopt) {
        auto tiledb_type = Object::object(*ctx->tiledb_ctx(), std::string(uri)).type();
        soma_type = util::soma_type_from_tiledb_type(tiledb_type);
    }

    if (soma_type == "SOMAArray") {
        auto array_ = SOMAArray::open(mode, uri, ctx, timestamp);
        auto array_type = array_->type();

        if (!array_type.has_value())
            throw TileDBSOMAError(fmt::format("[SOMAObject::open] SOMAArray '{}' has no type info", uri));

        std::transform(array_type->begin(), array_type->end(), array_type->begin(), [](unsigned char c) {
            return std::tolower(c);
        });

        if (array_type == "somadataframe") {
            return std::make_unique<SOMADataFrame>(*array_);
        } else if (array_type == "somasparsendarray") {
            return std::make_unique<SOMASparseNDArray>(*array_);
        } else if (array_type == "somadensendarray") {
            return std::make_unique<SOMADenseNDArray>(*array_);
        } else if (array_type == "somapointclouddataframe") {
            return std::make_unique<SOMAPointCloudDataFrame>(*array_);
        } else if (array_type == "somageometrydataframe") {
            return std::make_unique<SOMAGeometryDataFrame>(*array_);
        } else {
            throw TileDBSOMAError("[SOMAObject::open] Saw invalid SOMAArray type");
        }
    } else if (soma_type == "SOMAGroup") {
        auto group_ = SOMAGroup::open(mode, uri, ctx, "", timestamp);
        auto group_type = group_->type();

        if (!group_type.has_value())
            throw TileDBSOMAError(fmt::format("[SOMAObject::open] SOMAGroup '{}' has no type info", uri));

        std::transform(group_type->begin(), group_type->end(), group_type->begin(), [](unsigned char c) {
            return std::tolower(c);
        });

        if (group_type == "somacollection") {
            return std::make_unique<SOMACollection>(*group_);
        } else if (group_type == "somaexperiment") {
            return std::make_unique<SOMAExperiment>(*group_);
        } else if (group_type == "somameasurement") {
            return std::make_unique<SOMAMeasurement>(*group_);
        } else if (group_type == "somascene") {
            return std::make_unique<SOMAScene>(*group_);
        } else if (group_type == "somamultiscaleimage") {
            return std::make_unique<SOMAMultiscaleImage>(*group_);
        } else {
            throw TileDBSOMAError("[SOMAObject::open] Saw invalid SOMAGroup type");
        }
    }

    throw TileDBSOMAError("[SOMAObject::open] Invalid TileDB object passed to SOMAObject::open");
}

const std::optional<std::string> SOMAObject::type() {
    auto soma_object_type = this->get_metadata(SOMA_OBJECT_TYPE_KEY);

    if (!soma_object_type.has_value())
        return std::nullopt;

    const char* dtype = (const char*)std::get<MetadataInfo::value>(*soma_object_type);
    uint32_t sz = std::get<MetadataInfo::num>(*soma_object_type);

    return std::string(dtype, sz);
}

bool SOMAObject::check_type(std::string expected_type) {
    auto soma_object_type = this->type();

    if (!soma_object_type.has_value())
        return false;

    std::transform(soma_object_type->begin(), soma_object_type->end(), soma_object_type->begin(), [](unsigned char c) {
        return std::tolower(c);
    });

    std::transform(expected_type.begin(), expected_type.end(), expected_type.begin(), [](unsigned char c) {
        return std::tolower(c);
    });

    return soma_object_type == expected_type;
};

}  // namespace tiledbsoma
