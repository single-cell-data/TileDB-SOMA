/**
 * @file   common.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This declares the common functions in the API
 */

#ifndef TILEDBSOMA_COMMON_H
#define TILEDBSOMA_COMMON_H

#include <stdexcept>  // for windows: error C2039: 'runtime_error': is not a member of 'std'
#include <string>
#include <string_view>
#include <tiledb/tiledb>

namespace tiledbsoma {

constexpr std::string_view SOMA_JOINID = "soma_joinid";

const std::string SOMA_OBJECT_TYPE_KEY = "soma_object_type";
const std::string ENCODING_VERSION_KEY = "soma_encoding_version";
const std::string ENCODING_VERSION_VAL = "1.1.0";
const std::string
    SPATIAL_ENCODING_VERSION_KEY = "soma_spatial_encoding_version";
const std::string SPATIAL_ENCODING_VERSION_VAL = "0.2.0";
const std::string SOMA_COORDINATE_SPACE_KEY = "soma_coordinate_space";
const std::string SOMA_GEOMETRY_COLUMN_NAME = "soma_geometry";
const std::string SOMA_GEOMETRY_DIMENSION_PREFIX = "tiledb__internal__";
const std::string ARROW_DATATYPE_METADATA_KEY = "dtype";

// SOMAColumn metadata keys
const std::string TILEDB_SOMA_SCHEMA_KEY = "tiledb_soma_schema";
const std::string TILEDB_SOMA_SCHEMA_VERSION = "0.0.1";
const std::string TILEDB_SOMA_SCHEMA_COL_KEY = "tiledb_columns";
const std::string TILEDB_SOMA_SCHEMA_COL_TYPE_KEY = "tiledb_column_type";
const std::string TILEDB_SOMA_SCHEMA_COL_DIM_KEY = "tiledb_dimensions";
const std::string TILEDB_SOMA_SCHEMA_COL_ATTR_KEY = "tiledb_attributes";

using MetadataValue = std::tuple<tiledb_datatype_t, uint32_t, const void*>;
enum MetadataInfo { dtype = 0, num, value };

using TimestampRange = std::pair<uint64_t, uint64_t>;

class TileDBSOMAError : public std::runtime_error {
   public:
    explicit TileDBSOMAError(const char* m)
        : std::runtime_error(m){};
    explicit TileDBSOMAError(std::string m)
        : std::runtime_error(m.c_str()){};

   public:
    virtual const char* what() const noexcept override {
        return std::runtime_error::what();
    }
};

};  // namespace tiledbsoma

#endif  // TILEDBSOMA_COMMON_H
