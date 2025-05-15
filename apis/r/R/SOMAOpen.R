#' @title Open a SOMA Object
#' @description Utility function to open the corresponding SOMA Object given a URI, (lifecycle: maturing)
#' @param mode One of `"READ"` or `"WRITE"`
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @param tiledb_timestamp Optional Datetime (POSIXct) with TileDB timestamp. For SOMACollections,
#'        all accessed members inherit the collection opening timestamp, and in READ mode the
#'        collection timestamp defaults to the time of opening.
#'
#' @export
#'
SOMAOpen <- function(
  uri,
  mode = "READ",
  platform_config = NULL,
  tiledbsoma_ctx = NULL,
  tiledb_timestamp = NULL
) {
  ctx <- soma_context()
  metadata <- get_all_metadata(
    uri,
    is_array = switch(
      get_tiledb_object_type(uri, ctxxp = ctx),
      ARRAY = TRUE,
      GROUP = FALSE,
      stop("Unknown TileDB object type: ", dQuote(type), call. = FALSE)
    ),
    ctxxp = ctx
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
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMADataFrame = SOMADataFrameOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMADenseNDArray = SOMADenseNDArrayOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMASparseNDArray = SOMASparseNDArrayOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMAExperiment = SOMAExperimentOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    SOMAMeasurement = SOMAMeasurementOpen(
      uri,
      mode = mode,
      platform_config = platform_config,
      tiledbsoma_ctx = tiledbsoma_ctx,
      tiledb_timestamp = tiledb_timestamp
    ),
    stop(sprintf("No support for type '%s'", obj), call. = FALSE)
  ))
}
