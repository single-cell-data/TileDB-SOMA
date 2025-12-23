#' Create a SOMA Data Frame
#'
#' Factory function to create a \link[tiledbsoma:SOMADataFrame]{SOMA data frame}
#' for writing (lifecycle: maturing).
#'
#' @param uri URI for the TileDB object.
#' @param schema Arrow schema argument for the
#' \link[tiledbsoma:SOMADataFrame]{SOMA dataframe}.
#' @param index_column_names A vector of column names to use as user-defined
#' index columns; all named columns must exist in the schema, and at least
#' one index column name is required.
#' @param domain An optional list of 2-element vectors specifying the domain of
#' each index column. Each vector should be a pair consisting of the minimum and
#' maximum values storable in the index column. For example, if there is a
#' single int64-valued index column, then \code{domain} might be
#' \code{c(100, 200)} to indicate that values between 100 and 200, inclusive,
#' can be stored in that column. If provided, this list must have the same
#' length as \code{index_column_names}, and the index-column domain will be as
#' specified. If omitted entirely, or if \code{NULL} in a given dimension, the
#' corresponding index-column domain will use the minimum and maximum possible
#' values for the column's datatype. This makes a SOMA data frame growable.
#' @param ingest_mode Ingestion mode when creating the TileDB object;
#' choose from:
#' \itemize{
#'  \item \dQuote{\code{write}}: create a new TileDB object and error if
#'   it already exists.
#'  \item \dQuote{\code{resume}}: attempt to create a new TileDB object;
#'   if it already exists, simply open it for writing.
#' }
#' @param platform_config Optional platform configuration.
#' @param tiledbsoma_ctx Optional (DEPRECATED) SOMATileDBContext.
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp.
#' @param context Optional TileDB SOMA context.
#'
#' @return A new \link[tiledbsoma:SOMADataFrame]{SOMA data frame} stored at
#' \code{uri} opened for writing.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' uri <- withr::local_tempfile(pattern = "soma-data-frame")
#' df <- data.frame(
#'   soma_joinid = bit64::seq.integer64(0L, 99L),
#'   group = sample(factor(c("g1", "g2")), size = 100L, replace = TRUE),
#'   nCount = stats::rbinom(100L, 10L, 0.3)
#' )
#' (sch <- arrow::infer_schema(df))
#' (sdf <- SOMADataFrameCreate(uri, sch, domain = list(soma_joinid = c(0, 100))))
#' sdf$write(arrow::as_arrow_table(df, schema = sch))
#' sdf$close()
#'
#' (sdf <- SOMADataFrameOpen(uri))
#' head(as.data.frame(sdf$read()$concat()))
#' \dontshow{
#' sdf$close()
#' }
#'
SOMADataFrameCreate <- function(
  uri,
  schema,
  index_column_names = c("soma_joinid"),
  domain = NULL,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMADataFrameCreate(tiledbsoma_ctx)")
  sdf <- SOMADataFrame$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context,
    tiledb_timestamp = tiledb_timestamp
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(sdf$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    sdf$open(mode = "WRITE")
  } else {
    sdf$create(
      schema,
      index_column_names = index_column_names,
      domain = domain,
      platform_config = platform_config
    )
  }
  return(sdf)
}

#' Open a SOMA Data Frame
#'
#' Factory function to open a \link[tiledbsoma:SOMADataFrame]{SOMA data frame}
#' for reading (lifecycle: maturing).
#'
#' @inheritParams SOMADataFrameCreate
#' @param mode One of \dQuote{\code{READ}} or \dQuote{\code{WRITE}}.
#' @param platform_config Optional platform configuration.
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp;
#' defaults to the current time.
#'
#' @return A \link[tiledbsoma:SOMADataFrame]{SOMA data frame} stored at
#' \code{uri} opened in mode \code{mode}.
#'
#' @export
#'
#' @inherit SOMADataFrameCreate examples
#'
SOMADataFrameOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  soma_debug(sprintf(
    "[SOMADataFrameOpen] uri %s ts (%s)",
    uri,
    tiledb_timestamp %||% "now"
  ))
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMADataFrameOpen(tiledbsoma_ctx)")
  sdf <- SOMADataFrame$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context,
    tiledb_timestamp = tiledb_timestamp
  )
  sdf$open(mode)
  return(sdf)
}

#' Create a SOMA Sparse ND Array
#'
#' Factory function to create a
#' \link[tiledbsoma:SOMASparseNDArray]{SOMA sparse ND array} for writing
#' (lifecycle: maturing).
#'
#' @inheritParams SOMADataFrameCreate
#' @param type An \link[arrow:data-type]{Arrow type} defining the type
#' of each element in the array.
#' @param shape A vector of integers defining the shape of the array.
#'
#' @return A new \link[tiledbsoma:SOMASparseNDArray]{SOMA sparse ND array}
#' stored at \code{uri} opened for writing.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' uri <- withr::local_tempfile(pattern = "soma-sparse-array")
#' mat <- Matrix::rsparsematrix(100L, 100L, 0.7, repr = "T")
#' mat[1:3, 1:5]
#'
#' (arr <- SOMASparseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
#' arr$write(mat)
#' arr$close()
#'
#' (arr <- SOMASparseNDArrayOpen(uri))
#' m2 <- arr$read()$sparse_matrix()$concat()
#' m2[1:3, 1:5]
#' \dontshow{
#' arr$close()
#' }
#'
SOMASparseNDArrayCreate <- function(
  uri,
  type,
  shape,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMASparseNDArrayCreate(tiledbsoma_ctx)")
  snda <- SOMASparseNDArray$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    context = context,
    tiledb_timestamp = tiledb_timestamp
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(snda$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    snda$open(mode = "WRITE")
  } else {
    snda$create(type, shape, platform_config = platform_config)
  }
  return(snda)
}

#' Open a SOMA Sparse ND Array
#'
#' Factory function to open a
#' \link[tiledbsoma:SOMASparseNDArray]{SOMA sparse ND array} for reading
#' (lifecycle: maturing).
#'
#' @inheritParams SOMADataFrameOpen
#'
#' @return A \link[tiledbsoma:SOMASparseNDArray]{SOMA sparse ND array} stored at
#' \code{uri} opened in mode \code{mode}.
#'
#' @export
#'
#' @inherit SOMASparseNDArrayCreate examples
#'
SOMASparseNDArrayOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMASparseNDArrayOpen(tiledbsoma_ctx)")
  snda <- SOMASparseNDArray$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  snda$open(mode)
  return(snda)
}

#' Create a SOMA Dense ND Array
#'
#' Factory function to create a
#' \link[tiledbsoma:SOMADenseNDArray]{SOMA dense ND array} for writing
#' (lifecycle: maturing).
#'
#' @inheritParams SOMASparseNDArrayCreate
#'
#' @return A new \link[tiledbsoma:SOMADenseNDArray]{SOMA dense ND array}
#' stored at \code{uri} opened for writing.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' uri <- withr::local_tempfile(pattern = "soma-dense-array")
#' mat <- matrix(stats::rnorm(100L ^ 2L), nrow = 100L, ncol = 100L)
#' mat[1:3, 1:5]
#'
#' (arr <- SOMADenseNDArrayCreate(uri, arrow::float64(), shape = dim(mat)))
#' arr$write(mat)
#' arr$close()
#'
#' (arr <- SOMADenseNDArrayOpen(uri))
#' arr$read_arrow_table()
#' \dontshow{
#' arr$close()
#' }
#'
SOMADenseNDArrayCreate <- function(
  uri,
  type,
  shape,
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  soma_debug(sprintf(
    "[SOMADenseNDArrayCreate] tstamp (%s)",
    tiledb_timestamp %||% "now"
  ))
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMADenseNDArrayCreate(tiledbsoma_ctx)")
  dnda <- SOMADenseNDArray$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  dnda$create(type, shape, platform_config = platform_config)
  return(dnda)
}

#' Open a SOMA Dense Nd Array
#'
#' Factory function to open a
#' \link[tiledbsoma:SOMADenseNDArray]{SOMA dense ND array} for reading
#' (lifecycle: maturing).
#'
#' @inheritParams SOMADataFrameOpen
#'
#' @return A \link[tiledbsoma:SOMADenseNDArray]{SOMA dense ND array} stored at
#' \code{uri} opened in mode \code{mode}.
#'
#' @export
#'
#' @inherit SOMADenseNDArrayCreate examples
#'
SOMADenseNDArrayOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMADenseNDArrayOpen(tiledbsoma_ctx)")
  dnda <- SOMADenseNDArray$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  dnda$open(mode)
  return(dnda)
}

#' Create a SOMA Collection
#'
#' Factory function to create a
#' \link[tiledbsoma:SOMACollection]{SOMA collection} for writing
#' (lifecycle: maturing).
#'
#' @inheritParams SOMADataFrameCreate
#'
#' @return A new \link[tiledbsoma:SOMACollection]{SOMA collection} stored at
#' \code{uri} opened for writing.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' uri <- withr::local_tempfile(pattern = "soma-collection")
#'
#' (col <- SOMACollectionCreate(uri))
#' col$add_new_sparse_ndarray("sparse", arrow::float64(), shape = c(100L, 100L))
#' col$close()
#'
#' (col <- SOMACollectionOpen(uri))
#' col$names()
#'
#' \dontshow{
#' col$close()
#' }
#'
SOMACollectionCreate <- function(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMACollectionCreate(tiledbsoma_ctx)")
  coll <- SOMACollection$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(coll$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    coll$open(mode = "WRITE")
  } else {
    coll$create()
  }
  return(coll)
}

#' Open a SOMA Collection
#'
#' Factory function to open a \link[tiledbsoma:SOMACollection]{SOMA collection}
#' for reading (lifecycle: maturing).
#'
#' @inheritParams SOMADataFrameOpen
#' @param tiledb_timestamp Optional Datetime (POSIXct) for TileDB timestamp;
#' defaults to the current time. If not \code{NULL}, all members accessed
#' through the collection inherit the timestamp.
#'
#' @return A \link[tiledbsoma:SOMACollection]{SOMA collection} stored at
#' \code{uri} opened in mode \code{mode}.
#'
#' @export
#'
#' @inherit SOMACollectionCreate examples
#'
SOMACollectionOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMACollectionOpen(tiledbsoma_ctx)")
  coll <- SOMACollection$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  coll$open(mode)
  return(coll)
}

#' Create a SOMA Measurement
#'
#' Factory function to create a
#' \link[tiledbsoma:SOMAMeasurement]{SOMA measurement} for writing
#' (lifecycle: maturing).
#'
#' @inheritParams SOMACollectionCreate
#'
#' @return A new \link[tiledbsoma:SOMAMeasurement]{SOMA measurement} stored at
#' \code{uri} opened for writing.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' uri <- withr::local_tempfile(pattern = "soma-measurement")
#' var <- data.frame(
#'   soma_joinid = bit64::seq.integer64(0L, 99L),
#'   var_id = paste0("feature_", seq_len(100L))
#' )
#' sch <- arrow::infer_schema(var)
#'
#' (ms <- SOMAMeasurementCreate(uri))
#' sdf <- ms$add_new_dataframe(
#'   "var",
#'   sch,
#'   "soma_joinid",
#'   list(soma_joinid = c(0, 100))
#' )
#' sdf$write(arrow::as_arrow_table(var, schema = sch))
#' sdf$close()
#' ms$close()
#'
#' (ms <- SOMAMeasurementOpen(uri))
#' ms$var
#'
#' \dontshow{
#' ms$close()
#' }
#'
SOMAMeasurementCreate <- function(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMAMeasurementCreate(tiledbsoma_ctx)")
  meas <- SOMAMeasurement$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(meas$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    meas$open(mode = "WRITE")
  } else {
    meas$create()
  }
  return(meas)
}

#' Open SOMA Measurement
#'
#' Factory function to open a
#' \link[tiledbsoma:SOMAMeasurement]{SOMA measurement} for reading
#' (lifecycle: maturing).
#'
#' @inheritParams SOMACollectionOpen
#'
#' @return A \link[tiledbsoma:SOMAMeasurement]{SOMA measurement} stored at
#' \code{uri} opened in mode \code{mode}.
#'
#' @export
#'
#' @inherit SOMAMeasurementCreate examples
#'
SOMAMeasurementOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMAMeasurementOpen(tiledbsoma_ctx)")
  meas <- SOMAMeasurement$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  meas$open(mode)
  return(meas)
}

#' Create a SOMA Experiment
#'
#' Factory function to create a
#' \link[tiledbsoma:SOMAExperiment]{SOMA experiment} for writing
#' (lifecycle: maturing).
#'
#' @inheritParams SOMACollectionCreate
#'
#' @return A new \link[tiledbsoma:SOMAExperiment]{SOMA experiment} stored at
#' \code{uri} opened for writing.
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' uri <- withr::local_tempfile(pattern = "soma-experiment")
#' obs <- data.frame(
#'   soma_joinid = bit64::seq.integer64(0L, 99L),
#'   obs_id = paste0("cell_", seq_len(100L))
#' )
#' sch <- arrow::infer_schema(obs)
#'
#' (exp <- SOMAExperimentCreate(uri))
#' sdf <- exp$add_new_dataframe(
#'   "obs",
#'   sch,
#'   "soma_joinid",
#'   list(soma_joinid = c(0, 100))
#' )
#' sdf$write(arrow::as_arrow_table(obs, schema = sch))
#' sdf$close()
#' exp$close()
#'
#' (exp <- SOMAExperimentOpen(uri))
#' exp$obs
#'
#' \dontshow{
#' exp$close()
#' }
#'
SOMAExperimentCreate <- function(
  uri,
  ingest_mode = c("write", "resume"),
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  ingest_mode <- match.arg(ingest_mode)
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMAExperimentCreate(tiledbsoma_ctx)")
  exp <- SOMAExperiment$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  ingest_mode <- switch(
    EXPR = ingest_mode,
    resume = ifelse(exp$exists(), yes = ingest_mode, no = "write"),
    ingest_mode
  )
  if (ingest_mode %in% c("resume")) {
    exp$open(mode = "WRITE")
  } else {
    exp$create()
  }
  return(exp)
}

#' Open SOMA Experiment
#'
#' Factory function to open a \link[tiledbsoma:SOMAExperiment]{SOMA experiment}
#' for reading (lifecycle: maturing).
#'
#' @inheritParams SOMACollectionOpen
#'
#' @return A \link[tiledbsoma:SOMAExperiment]{SOMA experiment} stored at
#' \code{uri} opened in mode \code{mode}.
#'
#' @export
#'
#' @inherit SOMAExperimentCreate examples
#'
SOMAExperimentOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL,
  context = NULL
) {
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMAExperimentOpen(tiledbsoma_ctx)")
  exp <- SOMAExperiment$new(
    uri,
    platform_config = platform_config,
    tiledbsoma_ctx = tiledbsoma_ctx,
    tiledb_timestamp = tiledb_timestamp,
    context = context
  )
  exp$open(mode)
  return(exp)
}
