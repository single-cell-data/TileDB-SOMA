#' Open a SOMA Object
#'
#' Utility function to open the corresponding SOMA object given a URI
#' (lifecycle: maturing).
#'
#' @inheritParams SOMACollectionOpen
#' @param mode One of \dQuote{\code{READ}} or \dQuote{\code{WRITE}}
#'
#' @return A SOMA object
#'
#' @export
#'
#' @examplesIf requireNamespace("withr", quietly = TRUE)
#' dir <- withr::local_tempfile(pattern = "soma-open")
#' dir.create(dir, recursive = TRUE)
#'
#' uri <- extract_dataset("soma-exp-pbmc-small", dir)
#' (exp <- SOMAOpen(uri))
#'
#' \dontshow{
#' exp$close()
#' }
#'
#' uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs", dir)
#' (obs <- SOMAOpen(uri))
#'
#' \dontshow{
#' obs$close()
#' }
#'
SOMAOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  context = NULL,
  tiledb_timestamp = NULL
) {
  context <- get_soma_context(context, tiledbsoma_ctx, what="SOMAOpen(tiledbsoma_ctx)")
  type <- get_tiledb_object_type(uri, ctxxp = context$handle)
  metadata <- get_all_metadata(
    uri,
    is_array = switch(
      EXPR = type,
      ARRAY = TRUE,
      GROUP = FALSE,
      stop("Unknown TileDB object type: ", dQuote(type), call. = FALSE)
    ),
    ctxxp = context$handle
  )
  if (is.null(metadata$soma_object_type)) {
    stop("URI ", sQuote(uri), " is not a TileDB SOMA object", call. = FALSE)
  }

  return(switch(
    EXPR = metadata$soma_object_type,
    SOMACollection = SOMACollectionOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      context = context,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMADataFrame = SOMADataFrameOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      context = context,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMADenseNDArray = SOMADenseNDArrayOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      context = context,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMASparseNDArray = SOMASparseNDArrayOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      context = context,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMAExperiment = SOMAExperimentOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      context = context,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMAMeasurement = SOMAMeasurementOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      context = context,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    stop(
      "No support for type ",
      dQuote(metadata$soma_object_type),
      call. = FALSE
    )
  ))
}
