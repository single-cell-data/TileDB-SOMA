
#' @title Open a SOMA Object
#' @description Utility function to open the corresponding SOMA Object given a URI, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @export
SOMAOpen <- function(uri, platform_config = NULL, tiledbsoma_ctx = NULL) {
    # As an alternative we could rely tiledb-r and its tiledb_object_type but
    # this would require instantiating a ctx object first. It is a possible
    # refinement if and when we decide to hold an array or group pointer.
    # For now, first attempt to instantiate a TileDBArray to take advantage of
    # its handling of the config and ctx object as well as the caching
    obj <- tryCatch(expr = {
        arr <- TileDBArray$new(uri,
                               platform_config = platform_config,
                               tiledbsoma_ctx = tiledbsoma_ctx,
                               internal_use_only = "allowed_use")
        obj <- arr$get_metadata("soma_object_type")
    },
    error = function(...) NULL,
    finally = function(...) NULL)

    # In case of an error try again as TileDBGroup
    if (is.null(obj)) {
        obj <- tryCatch(expr = {
            arr <- TileDBGroup$new(uri,
                                   platform_config = platform_config,
                                   tiledbsoma_ctx = tiledbsoma_ctx,
                                   internal_use_only = "allowed_use")
            obj <- arr$get_metadata("soma_object_type")
        },
        error = function(...) NULL,
        finally = function(...) NULL)
    }

    # If this also errored no other
    if (is.null(obj)) {
        stop("URI '", uri, "' is not a TileDB SOMA object.", call. = FALSE)
    }

    # Alternative: go via tiledb-r but needs to set up platform_config + ctx first
    # TODO set to ctx from config and/or ctx for initial open too
    #arr <- tiledb::tiledb_array(uri) # sadly this currently returns it closed
    #arr <- tiledb::tiledb_array_open(arr, "READ") # TODO: just get it in opened state
    #obj <- tiledb::tiledb_get_metadata(arr, "soma_object_type")

    switch(obj,
           SOMACollection    = SOMACollectionOpen(uri, platform_config, tiledbsoma_ctx),
           SOMADataFrame     = SOMADataFrameOpen(uri, platform_config, tiledbsoma_ctx),
           SOMADenseNDArray  = SOMADenseNDArrayOpen(uri, platform_config, tiledbsoma_ctx),
           SOMASparseNDArray = SOMASparseNDArrayOpen(uri, platform_config, tiledbsoma_ctx),
           SOMAExperiment    = SOMAExperimentOpen(uri, platform_config, tiledbsoma_ctx),
           SOMAMeasurement   = SOMAMeasurementOpen(uri, platform_config, tiledbsoma_ctx),
           stop(sprintf("No support for type '%s'", obj), call. = FALSE))
}
