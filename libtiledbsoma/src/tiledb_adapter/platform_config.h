/**
 * @file   platform_config.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 * This file defines classes and helper methods for the TileDB platform configuration.
 */

#ifndef SOMA_PLATFORM_CONFIG_H
#define SOMA_PLATFORM_CONFIG_H

#include <cstdint>
#include <optional>
#include <span>
#include <string>

#include <tiledb/tiledb>

#include "nlohmann/json.hpp"

namespace tiledbsoma {

using namespace tiledb;

struct PlatformConfig {
   public:
    /* Set the ZstdFilter's level for DataFrame dims */
    int32_t dataframe_dim_zstd_level = 3;

    /* Set the ZstdFilter's level for SparseNDArray dims */
    int32_t sparse_nd_array_dim_zstd_level = 3;

    /* Set the ZstdFilter's level for DenseNDArray dims */
    int32_t dense_nd_array_dim_zstd_level = 3;

    /* Set whether to write the X data chunked */
    bool write_X_chunked = true;

    /* Set the goal chunk nnz */
    uint64_t goal_chunk_nnz = 100000000;

    /* Server-side parameter to set the cap nbytes */
    uint64_t remote_cap_nbytes = 2400000000;

    /* Set the tile capcity for sparse arrays */
    uint64_t capacity = 100000;

    /**
     * Available filters with associated options are
     * [
     *     {
     *         "name": "GZIP", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "ZSTD", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "LZ4", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "BZIP2", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "RLE", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "DELTA",
     *         "COMPRESSION_LEVEL": (int32_t),
     *         "COMPRESSION_REINTERPRET_DATATYPE": (uint8_t)
     *     },
     *     {
     *         "name": "DOUBLE_DELTA",
     *         "COMPRESSION_LEVEL": (int32_t),
     *         "COMPRESSION_REINTERPRET_DATATYPE": (uint8_t)
     *     },
     *     {
     *         "name": "BIT_WIDTH_REDUCTION",
     *         "BIT_WIDTH_MAX_WINDOW": (uint32_t)
     *     },
     *     {
     *         "name": "POSITIVE_DELTA", "POSITIVE_DELTA_MAX_WINDOW":
     * (uint32_t),
     *     },
     *     {
     *         "name": "DICTIONARY_ENCODING", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "SCALE_FLOAT",
     *         "SCALE_FLOAT_FACTOR": (double),
     *         "SCALE_FLOAT_OFFSET": (double),
     *         "SCALE_FLOAT_BYTEWIDTH": (uint64_t),
     *     },
     *     {
     *         "name": "WEBP",
     *         "WEBP_INPUT_FORMAT": (uint8_t),
     *         "WEBP_QUALITY": (float),
     *         "WEBP_LOSSLESS": (uint8_t),
     *     },
     *     "CHECKSUM_MD5",
     *     "CHECKSUM_SHA256",
     *     "XOR",
     *     "BITSHUFFLE",
     *     "BYTESHUFFLE",
     *     "NOOP"
     * ]
     *
     */
    std::string offsets_filters = R"(["DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD"])";

    /* Set the validity filters. */
    std::string validity_filters = "";

    /* Set whether the TileDB Array allows duplicate values */
    bool allows_duplicates = false;

    /* Set the tile order as "row", "row-major", "col", or "col-major" */
    std::optional<std::string> tile_order = std::nullopt;

    /* Set the cell order as "hilbert", "row", "row-major", "col", or
     * "col-major"
     */
    std::optional<std::string> cell_order = std::nullopt;

    /* Set the filters for attributes.
     *
     * Example:
     * {
     *     "attr_name": {
     *          "filters": ["XOR", {"name": "GZIP", "COMPRESSION_LEVEL": 3}]
     *     }
     * }
     *
     */
    std::string attrs = "";

    /* Set the filters and tiles for dimensions.
     *
     * Example:
     * {
     *     "dim_name": {"filters": ["NoOpFilter"], "tile": 8}
     * }
     *
     */

    std::string dims = "";

    /* Set whether the array should be consolidated and vacuumed after writing
     */
    bool consolidate_and_vacuum = false;
};

/** TileDB specific configuration options that can be read back from a single
 * TileDB ArraySchema.
 */
struct PlatformSchemaConfig {
   public:
    /* Set whether the TileDB Array allows duplicate values */
    bool allows_duplicates = false;

    /* Set the tile order as "row", "row-major", "col", or "col-major" */
    std::optional<std::string> tile_order = std::nullopt;

    /* Set the cell order as "hilbert", "row", "row-major", "col", or
     * "col-major"
     */
    std::optional<std::string> cell_order = std::nullopt;

    /* Set the tile capcity for sparse arrays */
    uint64_t capacity = 100000;

    /**
     * Available filters with associated options are
     * [
     *     {
     *         "name": "GZIP", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "ZSTD", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "LZ4", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "BZIP2", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "RLE", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "DELTA",
     *         "COMPRESSION_LEVEL": (int32_t),
     *         "COMPRESSION_REINTERPRET_DATATYPE": (uint8_t)
     *     },
     *     {
     *         "name": "DOUBLE_DELTA",
     *         "COMPRESSION_LEVEL": (int32_t),
     *         "COMPRESSION_REINTERPRET_DATATYPE": (uint8_t)
     *     },
     *     {
     *         "name": "BIT_WIDTH_REDUCTION",
     *         "BIT_WIDTH_MAX_WINDOW": (uint32_t)
     *     },
     *     {
     *         "name": "POSITIVE_DELTA", "POSITIVE_DELTA_MAX_WINDOW":
     * (uint32_t),
     *     },
     *     {
     *         "name": "DICTIONARY_ENCODING", "COMPRESSION_LEVEL": (int32_t)
     *     },
     *     {
     *         "name": "SCALE_FLOAT",
     *         "SCALE_FLOAT_FACTOR": (double),
     *         "SCALE_FLOAT_OFFSET": (double),
     *         "SCALE_FLOAT_BYTEWIDTH": (uint64_t),
     *     },
     *     {
     *         "name": "WEBP",
     *         "WEBP_INPUT_FORMAT": (uint8_t),
     *         "WEBP_QUALITY": (float),
     *         "WEBP_LOSSLESS": (uint8_t),
     *     },
     *     "CHECKSUM_MD5",
     *     "CHECKSUM_SHA256",
     *     "XOR",
     *     "BITSHUFFLE",
     *     "BYTESHUFFLE",
     *     "NOOP"
     * ]
     *
     */
    std::string offsets_filters = R"(["DOUBLE_DELTA", "BIT_WIDTH_REDUCTION", "ZSTD"])";

    /* Set the validity filters. */
    std::string validity_filters = "";

    /* Set the filters for attributes.
     *
     * Example:
     * {
     *     "attr_name": {
     *          "filters": ["XOR", {"name": "GZIP", "COMPRESSION_LEVEL": 3}]
     *     }
     * }
     *
     */
    std::string attrs = "";

    /* Set the filters and tiles for dimensions.
     *
     * Example:
     * {
     *     "dim_name": {"filters": ["NoOpFilter"], "tile": 8}
     * }
     *
     */

    std::string dims = "";
};

template <typename T>
class DimensionConfigAdapter {
   public:
    /**
    * Constructor using Arrow buffer from API.
    *
    * Buffer store the following values.
    * 0: max domain lower bound
    * 1: max domain upper bound
    * 2: tile extent
    * 3: current domain lower bound
    * 4: current domain upper bound
    *
    * @parm buffer Buffer storing configuration paramters.
    *
    */
    DimensionConfigAdapter(const T* buffer)
        : current_domain_{buffer[3], buffer[4]}
        , max_domain_{buffer[0], buffer[1]}
        , tile_extent_{buffer[2]} {
    }

    DimensionConfigAdapter(std::array<T, 2> current_domain, std::array<T, 2> max_domain, T tile_extent)
        : current_domain_{current_domain}
        , max_domain_{max_domain}
        , tile_extent_{tile_extent} {
    }

    std::array<T, 2> max_domain() const {
        return max_domain_;
    }

    std::array<T, 2> current_domain() const {
    }

    Dimension create_dimension(const Context& ctx, std::string name, const FilterList& filter_list) const {
        auto dim = Dimension::create(ctx, name, max_domain_, tile_extent_);
        dim.set_filter_list(filter_list);
        return dim;
    }

    Dimension create_dimension(
        const Context& ctx, std::string name, tiledb_datatype_t datatype, const FilterList& filter_list) const {
        auto dim = Dimension::create(ctx, name, datatype, max_domain_.data(), &tile_extent_);
        dim.set_filter_list(filter_list);
        return dim;
    }

   private:
    std::array<T, 2> current_domain_;
    std::array<T, 2> max_domain_;
    T tile_extent_;
};

namespace utils {

ArraySchema create_base_tiledb_schema(
    std::shared_ptr<Context> ctx,
    const PlatformConfig& platform_config,
    bool is_sparse,
    std::optional<std::pair<int64_t, int64_t>> timestamp_range = std::nullopt);

/**
 * Create a TileDB filter list from JSON serialized in a string.
 *
 * @param filters A json description of filters serialized to a string.
 * @param ctx TileDB context.
 * @return FilterList
 */
FilterList create_filter_list(std::string filters, std::shared_ptr<Context> ctx);

/**
 * Create a TileDB filter list from JSON serialized in a string.
 *
 * @param name The name of the TileDB attribute.
 * @param platform_config A platform config object for the array the will contain the attribute.
 * @param ctx TileDB context.
 * @return FilterList
 */
FilterList create_attr_filter_list(std::string name, PlatformConfig platform_config, std::shared_ptr<Context> ctx);

/**
 * Create a TileDB filter list from JSON serialized in a string.
 *
 * @param name The name of the TileDB dimension.
 * @param soma_type The type of the SOMA object that will be created.
 * @param platform_config A platform config object for the array the will contain the attribute.
 * @param ctx TileDB context.
 */
FilterList create_dim_filter_list(
    std::string name, PlatformConfig platform_config, std::string soma_type, std::shared_ptr<Context> ctx);

template <typename T>
T get_dim_extent(std::string name, PlatformConfig platform_config, T default_extent, T max_domain) {
    T extent;
    if (platform_config.dims.empty()) {
        extent = default_extent;
    } else {
        nlohmann::json dim_options = nlohmann::json::parse(platform_config.dims);
        if (dim_options.find(name) != dim_options.end() && dim_options[name].find("tile") != dim_options[name].end()) {
            extent = dim_options[name].value("tile", default_extent);
        } else {
            extent = default_extent;
        }
    }

    return std::min(extent, max_domain);
}

/**
* Get members of the TileDB Schema in the form of a PlatformSchemaConfig.
*
* @param tiledb_schema The TileDB Schema to convert to a PlatformSchemaConfig.
* @return PlatformSchemaConfig
*/
PlatformSchemaConfig platform_schema_config_from_tiledb(ArraySchema tiledb_schema);

/**
* Get members of the TileDB Schema in the form of a PlatformConfig.
*
* @param tiledb_schema The TileDB Schema to convert to a PlatformConfig.
* @return PlatformConfig
*/
PlatformConfig platform_config_from_tiledb_schema(ArraySchema tiledb_schema);

ArraySchema create_nd_array_schema(
    std::string_view soma_type,
    bool is_sparse,
    std::string_view format,
    std::span<const int64_t> shape,
    std::shared_ptr<tiledb::Context> ctx,
    PlatformConfig platform_config,
    std::optional<std::pair<uint64_t, uint64_t>> timestamp);

}  // namespace utils

}  // namespace tiledbsoma

#endif
