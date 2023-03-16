
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
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @export
SOMADataFrameOpen <- function(uri, platform_config = NULL, tiledbsoma_ctx = NULL) {
    sdf <- SOMADataFrame$new(uri, platform_config, tiledbsoma_ctx, mode="READ",
                             internal_use_only = "allowed_use")
    ## TODO: other things to cache ?
    ## TODO: explicitly open and hold handle ?
    sdf
}

#' @title Create SOMA Sparse Nd Array
#' @description Factory function to create a SOMASparseNDArray for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB objec
#' @param type An [Arrow type][arrow::data-type] defining the type of each element in the array.
#' @param shape A vector of integers defining the shape of the array.
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
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @export
SOMASparseNDArrayOpen <- function(uri, platform_config = NULL, tiledbsoma_ctx = NULL) {
    sdar <- SOMASparseNDArray$new(uri, platform_config, tiledbsoma_ctx,
                                  mode="READ", internal_use_only = "allowed_use")
    ## TODO: other things to cache ?
    ## TODO: explicitly open and hold handle ?
    sdar
}

#' @title Create SOMA Dense Nd Array
#' @description Factory function to create a SOMADenseNDArray for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB objec
#' @param type An [Arrow type][arrow::data-type] defining the type of each element in the array.
#' @param shape A vector of integers defining the shape of the array.
#' @export
SOMADenseNDArrayCreate <- function(uri, type, shape) {
    dnar <- SOMADenseNDArray$new(uri, mode="WRITE", internal_use_only = "allowed_use")
    dnar$create(type, shape)
    dnar
}

#' @title Open SOMA Dense Nd Array
#' @description Factory function to open a SOMADenseNDArray for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @export
SOMADenseNDArrayOpen <- function(uri, platform_config = NULL, tiledbsoma_ctx = NULL) {
    dnar <- SOMADenseNDArray$new(uri, platform_config, tiledbsoma_ctx,
                                 mode="READ", internal_use_only = "allowed_use")
    ## TODO: other things to cache ?
    ## TODO: explicitly open and hold handle ?
    dnar
}

#' @title Create SOMA Collection
#' @description Factory function to create a SOMADataFrame for writing, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx Optional SOMATileDBContext
#' @export
SOMACollectionCreate <- function(uri, platform_config = NULL, tiledbsoma_ctx = NULL) {
    clct <- SOMACollection$new(uri, platform_config, tiledbsoma_ctx,
                               internal_use_only = "allowed_use")
    clct$create()
    clct
}

#' @title Open SOMA Collection
#' @description Factory function to open a SOMACollection for reading, (lifecycle: experimental)
#' @param uri URI for the TileDB object
#' @param platform_config Optional platform configuration
#' @param tiledbsoma_ctx optional SOMATileDBContext
#' @export
SOMACollectionOpen <- function(uri, platform_config = NULL, tiledbsoma_ctx = NULL) {
    sdf <- SOMACollection$new(uri, platform_config, tiledbsoma_ctx,
                             internal_use_only = "allowed_use")
    ## TODO: other things to cache ?
    ## TODO: explicitly open and hold handle ?
    sdf
}
