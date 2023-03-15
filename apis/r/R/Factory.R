
#' @title Create SOMA DataFrame
#' @description Factory function to create a SOMADataFrame for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param schema schema Arrow schema argument passed on to DataFrame$create()
#' @param index_column_names Index column names passed on to DataFrame$create()
#' @export
SOMADataFrameCreate <- function(uri, schema, index_column_names) {
    sdf <- SOMADataFrame$new(uri, mode="WRITE", internal_use_only = "allowed_use")
    sdf$create(schema, index_column_names)
    sdf
}

#' @title Open SOMA DataFrame
#' @description Factory function to open a SOMADataFrame for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx optional SOMATileDBContext
#' @export
SOMADataFrameOpen <- function(uri, platform_config = NULL, tiledbsoma_ctx = NULL) {
    sdf <- SOMADataFrame$new(uri, platform_config, tiledbsoma_ctx, mode="READ",
                             internal_use_only = "allowed_use")
    ## TODO: other things to cache ?
    ## TODO: explicitly open and hold handle ?
    sdf
}

#' @title Create SOMA Sparse Nd Array
#' @description Factory function to create a SOMASparseNdArray for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB objec
#' @param type an [Arrow type][arrow::data-type] defining the type of each element in the array.
#' @param shape a vector of integers defining the shape of the array.
#' @export
SOMASparseNDArrayCreate <- function(uri, type, shape) {
    spar <- SOMASparseNDArray$new(uri, mode="WRITE", internal_use_only = "allowed_use")
    spar$create(type, shape)
    spar
}

#' @title Open SOMA Sparse Nd Array
#' @description Factory function to open a SOMASparseNDArray for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx optional SOMATileDBContext
#' @export
SOMASparseNDArrayOpen <- function(uri, platform_config = NULL, tiledbsoma_ctx = NULL) {
    sdf <- SOMASparseNDArray$new(uri, platform_config, tiledbsoma_ctx,
                                 mode="READ", internal_use_only = "allowed_use")
    ## TODO: other things to cache ?
    ## TODO: explicitly open and hold handle ?
    sdf
}
