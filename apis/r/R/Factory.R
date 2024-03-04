
#' @title Create SOMA DataFrame
#' @description Factory function to create a SOMADataFrame for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param schema schema Arrow schema argument passed on to DataFrame$create()
#' @param index_column_names Index column names passed on to DataFrame$create()
#' @param mode Ingestion mode: one of \code{write} or \code{resume}
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMADataFrameCreate <- function(uri, schema, index_column_names = c("soma_joinid"),
                                mode = "write",
                                platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    spdl::debug("[SOMADataFrameCreate] mode={}", mode)
    sdf <- SOMADataFrame$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp, internal_use_only = "allowed_use")
    if (!(mode == "resume" && tiledb::tiledb_vfs_is_dir(sdf$uri))) {
      sdf$create(schema, index_column_names=index_column_names,
                 platform_config=platform_config, internal_use_only = "allowed_use")
    }
    sdf
}

#' @title Open SOMA DataFrame
#' @description Factory function to open a SOMADataFrame for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param mode One of `"READ"` or `"WRITE"`
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMADataFrameOpen <- function(uri, mode="READ",
                              platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    sdf <- SOMADataFrame$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                             internal_use_only = "allowed_use")
    sdf$open(mode, internal_use_only = "allowed_use")
    sdf
}

#' @title Create SOMA Sparse Nd Array
#' @description Factory function to create a SOMASparseNDArray for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param type An [Arrow type][arrow::data-type] defining the type of each element in the array.
#' @param shape A vector of integers defining the shape of the array.
#' @param mode Ingestion mode: one of \code{write} or \code{resume}
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMASparseNDArrayCreate <- function(uri, type, shape, mode = "write",
                                    platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    spdl::debug("[SOMASparseArrayCreate] mode={}", mode)
    snda <- SOMASparseNDArray$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                                  internal_use_only = "allowed_use")
    spdl::debug("[SOMASparseArrayCreate] array at {}", snda$uri)
    if (!(mode == "resume" && tiledb::tiledb_vfs_is_dir(snda$uri))) {
        spdl::debug("[SOMASparseArrayCreate] creating uri={}", snda$uri)
        snda$create(type, shape, platform_config=platform_config, internal_use_only = "allowed_use")
    }
    snda
}

#' @title Open SOMA Sparse Nd Array
#' @description Factory function to open a SOMASparseNDArray for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param mode One of `"READ"` or `"WRITE"`
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMASparseNDArrayOpen <- function(uri, mode="READ",
                                  platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    snda <- SOMASparseNDArray$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                                  internal_use_only = "allowed_use")
    snda$open(mode, internal_use_only = "allowed_use")
    snda
}

#' @title Create SOMA Dense Nd Array
#' @description Factory function to create a SOMADenseNDArray for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param type An [Arrow type][arrow::data-type] defining the type of each element in the array.
#' @param shape A vector of integers defining the shape of the array.
#' @param mode Ingestion mode: one of \code{write} or \code{resume}
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMADenseNDArrayCreate <- function(uri, type, shape, mode = "write",
                                   platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    spdl::debug("[SOMADenseArrayCreate] mode={}", mode)
    dnda <- SOMADenseNDArray$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                                 internal_use_only = "allowed_use")
    if (!(mode == "resume" && tiledb::tiledb_vfs_is_dir(dnda$uri))) {
        spdl::debug("[SOMADenseArrayCreate] creating uri={}", dnda$uri)
        dnda$create(type, shape, platform_config=platform_config, internal_use_only = "allowed_use")
    }
    dnda
}

#' @title Open SOMA Dense Nd Array
#' @description Factory function to open a SOMADenseNDArray for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param mode One of `"READ"` or `"WRITE"`
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMADenseNDArrayOpen <- function(uri, mode="READ",
                                 platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    dnda <- SOMADenseNDArray$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                                 internal_use_only = "allowed_use")
    dnda$open(mode, internal_use_only = "allowed_use")
    dnda
}

#' @title Create SOMA Collection
#' @description Factory function to create a SOMADataFrame for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMACollectionCreate <- function(uri,
                                 platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    coll <- SOMACollection$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                               internal_use_only = "allowed_use")
    coll$create(internal_use_only = "allowed_use")
    coll
}

#' @title Open SOMA Collection
#' @description Factory function to open a SOMACollection for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param mode One of `"READ"` or `"WRITE"`
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp. In READ mode, defaults
#'        to the current time. If non-NULL, then all members accessed through the collection object
#'        inherit the timestamp.
#' @export
SOMACollectionOpen <- function(uri, mode="READ",
                               platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    coll <- SOMACollection$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                               internal_use_only = "allowed_use")
    coll$open(mode, internal_use_only = "allowed_use")
    coll
}

#' @title Create SOMA Measurement
#' @description Factory function to create a SOMAMeasurement for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMAMeasurementCreate <- function(uri,
                                  platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    meas <- SOMAMeasurement$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                                internal_use_only = "allowed_use")
    meas$create(internal_use_only = "allowed_use")
    meas
}

#' @title Open SOMA Measurement
#' @description Factory function to open a SOMAMeasurement for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param mode One of `"READ"` or `"WRITE"`
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp. In READ mode, defaults
#'        to the current time. If non-NULL, then all members accessed through the collection object
#'        inherit the timestamp.
#' @export
SOMAMeasurementOpen <- function(uri, mode="READ",
                                platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    meas <- SOMAMeasurement$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                                internal_use_only = "allowed_use")
    meas$open(mode, internal_use_only = "allowed_use")
    meas
}

#' @title Create SOMA Experiment
#' @description Factory function to create a SOMADataFrame for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#' @export
SOMAExperimentCreate <- function(uri,
                                 platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    exp <- SOMAExperiment$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                              internal_use_only = "allowed_use")
    exp$create(internal_use_only = "allowed_use")
    exp
}

#' @title Open SOMA Experiment
#' @description Factory function to open a SOMAExperiment for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param mode One of `"READ"` or `"WRITE"`
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp. In READ mode, defaults
#'        to the current time. If non-NULL, then all members accessed through the collection object
#'        inherit the timestamp.
#' @export
SOMAExperimentOpen <- function(uri, mode="READ",
                               platform_config = NULL, tiledbsoma_ctx = NULL, tiledb_timestamp = NULL) {
    exp <- SOMAExperiment$new(uri, platform_config, tiledbsoma_ctx, tiledb_timestamp,
                              internal_use_only = "allowed_use")
    exp$open(mode, internal_use_only = "allowed_use")
    exp
}
