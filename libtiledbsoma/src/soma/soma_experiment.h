/**
 * @file   soma_experiment.h
 *
 * @section LICENSE
 *
 * Licensed under the MIT License.
 * Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
 *
 * @section DESCRIPTION
 *
 *   This file defines the SOMAExperiment class.
 */

#ifndef SOMA_EXPERIMENT
#define SOMA_EXPERIMENT

#include <tiledb/tiledb>

#include "soma_collection.h"
#include "soma_dataframe.h"

namespace tiledbsoma {

using namespace tiledb;
class SOMAExperiment : public SOMACollection {
   public:
    //===================================================================
    //= public static
    //===================================================================

    /**
     * @brief Create a SOMAExperiment object at the given URI.
     *
     * @param uri URI to create the SOMAExperiment
     * @param schema TileDB ArraySchema
     * @param platform_config Optional config parameter dictionary
     */
    static void create(
        std::string_view uri,
        const managed_unique_ptr<ArrowSchema>& schema,
        const ArrowTable& index_columns,
        std::shared_ptr<SOMAContext> ctx,
        PlatformConfig platform_config = PlatformConfig(),
        std::optional<TimestampRange> timestamp = std::nullopt);

    /**
     * @brief Open a group at the specified URI and return SOMAExperiment
     * object.
     *
     * @param uri URI of the array
     * @param mode read or write
     * @param ctx TileDB context
     * @param timestamp Optional pair indicating timestamp start and end
     * @return std::shared_ptr<SOMAExperiment> SOMAExperiment
     */
    static std::unique_ptr<SOMAExperiment> open(
        std::string_view uri,
        OpenMode mode,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt);

    //===================================================================
    //= public non-static
    //===================================================================

    SOMAExperiment(
        OpenMode mode,
        std::string_view uri,
        std::shared_ptr<SOMAContext> ctx,
        std::optional<TimestampRange> timestamp = std::nullopt)
        : SOMACollection(mode, uri, ctx, timestamp) {
    }

    SOMAExperiment(const SOMACollection& other)
        : SOMACollection(other) {
    }

    SOMAExperiment() = delete;
    SOMAExperiment(const SOMAExperiment&) = default;
    SOMAExperiment(SOMAExperiment&&) = default;
    ~SOMAExperiment() = default;

    /**
     * @brief Get the primary annotations on the observation axis
     * @param column_names A list of column names to use as user-defined
     index
     * columns (e.g., ``['cell_type', 'tissue_type']``). All named columns
     must
     * exist in the schema, and at least one index column name is required.
     * @param result_order Read result order: automatic (default), rowmajor,
     or
     * colmajor
     *
     * @return std::shared_ptr<SOMADataFrame>
     */
    std::shared_ptr<SOMADataFrame> obs();

    /**
     * @brief Get the collection of named measurements
     *
     * @return std::shared_ptr<SOMACollection>
     */
    std::shared_ptr<SOMACollection> ms();

   private:
    //===================================================================
    //= private non-static
    //===================================================================

    // Primary annotations on the observation axis
    std::shared_ptr<SOMADataFrame> obs_ = nullptr;

    // A collection of named measurements
    std::shared_ptr<SOMACollection> ms_ = nullptr;

    // A collection of spatial scenes
    std::shared_ptr<SOMACollection> spatial_ = nullptr;
};
}  // namespace tiledbsoma

#endif  // SOMA_EXPERIMENT
