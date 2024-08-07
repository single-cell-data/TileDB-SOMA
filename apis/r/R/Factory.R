
#' Create SOMA DataFrame
#'
#' Factory function to create a SOMADataFrame for writing, (lifecycle: maturing)
#'
#' @param uri URI for the TileDB object
#' @param schema Arrow schema argument for the \link[SOMADataFrame]{SOMA dataframe}
#' @param index_column_names A vector of column names to use as user-defined
#' index columns; all named columns must exist in the schema, and at least
#' one index column name is required
#' @param ingest_mode Ingestion mode when creating the TileDB object; choose from:
#' \itemize{
#'  \item \dQuote{\code{write}}: create a new TileDB object and error if it already exists
#'  \item \dQuote{\code{resume}}: attempt to create a new TileDB object;
#'   if it already exists, simply open it for writing
#' }
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp
#'
#' @export
#'
SOMADataFrameCreate <- function(
  uri,
  schema,
  index_column_names = c("soma_joinid"),
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  sdf <- SOMADataFrame$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(sdf$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    sdf$open(mode = "WRITE", internal_use_only = "allowed_use")
  } else {
    sdf$create(
      schema,
      index_column_names = index_column_names,
      platform_config = platform_config,
      timestamps = rep(tiledb_timestamp, 2),
      internal_use_only = "allowed_use"
    )
  }
  return(sdf)
}

#' Open SOMA DataFrame
#'
#' Factory function to open a SOMADataFrame for reading, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameCreate
#' @param mode One of `"READ"` or `"WRITE"`
#' @param platform_config Optional platform configuration
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp.
#' In READ mode, defaults to the current time. If non-NULL, then all members
#' accessed through the collection object inherit the timestamp.
#'
#' @export
#'
SOMADataFrameOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  sdf <- SOMADataFrame$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  sdf$open(mode, internal_use_only = "allowed_use")
  return(sdf)
}

#' Create SOMA Sparse Nd Array
#'
#' Factory function to create a SOMASparseNDArray for writing, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameCreate
#' @param type An [Arrow type][arrow::data-type] defining the type of each element in the array.
#' @param shape A vector of integers defining the shape of the array.
#'
#' @export
#'
SOMASparseNDArrayCreate <- function(
  uri,
  type,
  shape,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  snda <- SOMASparseNDArray$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(snda$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    snda$open(mode = "WRITE", internal_use_only = "allowed_use")
  } else {
    snda$create(
      type,
      shape,
      platform_config = platform_config,
      timestamps = rep(tiledb_timestamp, 2),
      internal_use_only = "allowed_use"
    )
  }
  return(snda)
}

#' Open SOMA Sparse Nd Array
#'
#' Factory function to open a SOMASparseNDArray for reading, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameOpen
#'
#' @export
#'
SOMASparseNDArrayOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  snda <- SOMASparseNDArray$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  snda$open(mode, internal_use_only = "allowed_use")
  return(snda)
}

#' Create SOMA Dense Nd Array
#'
#' Factory function to create a SOMADenseNDArray for writing, (lifecycle: maturing)
#'
#' @inheritParams SOMASparseNDArrayCreate
#'
#' @export
#'
SOMADenseNDArrayCreate <- function(
  uri,
  type,
  shape,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  dnda <- SOMADenseNDArray$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  dnda$create(
    type,
    shape,
    platform_config = platform_config,
    timestamps = rep(tiledb_timestamp, 2),
    internal_use_only = "allowed_use"
  )
  return(dnda)
}

#' Open SOMA Dense Nd Array
#'
#' Factory function to open a SOMADenseNDArray for reading, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameOpen
#'
#' @export
#'
SOMADenseNDArrayOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  dnda <- SOMADenseNDArray$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  dnda$open(mode, internal_use_only = "allowed_use")
  return(dnda)
}

#' Create SOMA Collection
#'
#' Factory function to create a SOMADataFrame for writing, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameCreate
#'
#' @export
#'
SOMACollectionCreate <- function(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  coll <- SOMACollection$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(coll$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    coll$open(mode = "WRITE", internal_use_only = "allowed_use")
  } else {
    coll$create(internal_use_only = "allowed_use")
  }
  return(coll)
}

#' Open SOMA Collection
#'
#' Factory function to open a SOMACollection for reading, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameOpen
#'
#' @export
#'
SOMACollectionOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  coll <- SOMACollection$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  coll$open(mode, internal_use_only = "allowed_use")
  return(coll)
}

#' Create SOMA Measurement
#'
#' Factory function to create a SOMAMeasurement for writing, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameCreate
#'
#' @export
#'
SOMAMeasurementCreate <- function(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  meas <- SOMAMeasurement$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(meas$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    meas$open(mode = "WRITE", internal_use_only = "allowed_use")
  } else {
    meas$create(internal_use_only = "allowed_use")
  }
  return(meas)
}

#' Open SOMA Measurement
#'
#' Factory function to open a SOMAMeasurement for reading, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameOpen
#'
#' @export
#'
SOMAMeasurementOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  meas <- SOMAMeasurement$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  meas$open(mode, internal_use_only = "allowed_use")
  return(meas)
}

#' Create SOMA Experiment
#'
#' Factory function to create a SOMADataFrame for writing, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameCreate
#'
#' @export
#'
SOMAExperimentCreate <- function(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  exp <- SOMAExperiment$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(exp$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    exp$open(mode = "WRITE", internal_use_only = "allowed_use")
  } else {
    exp$create(internal_use_only = "allowed_use")
  }
  return(exp)
}

#' Open SOMA Experiment
#'
#' Factory function to open a SOMAExperiment for reading, (lifecycle: maturing)
#'
#' @inheritParams SOMADataFrameOpen
#'
#' @export
#'
SOMAExperimentOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  exp <- SOMAExperiment$new(
    uri,
    platform_config,
    tiledbsoma_ctx,
    tiledb_timestamp,
    internal_use_only = "allowed_use"
  )
  exp$open(mode, internal_use_only = "allowed_use")
  return(exp)
}
