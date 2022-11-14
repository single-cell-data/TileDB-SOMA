##' Prototype Interfaces to \code{libtiledbsoma}
##'
##' Access functionality of \code{libtiledbsoma}
##'
##' @param uri A character variable with the URI of a TileDB (SOMA) Array
##' @param column,colname,colnames A character variable with a column name
##' @param arr A \code{arch} \code{array} with \code{arrow} representation of a \code{RecordBatch}
##' @param schemaxp A external pointer to an \code{arrow} schema object
##' @param arrayxp A external pointer to an \code{arrow} array data object
##' @param level A character variable with the desired logging level
##' @return The different functions return their respective values
get_table <- function(uri) {
    colnames <- get_column_names(uri)
    ll <- lapply(colnames, function(n) get_column(uri, n))
    names(ll) <- colnames
    as.data.frame(ll)
}

##' @rdname get_table
get_column <- function(uri, column) {
    schema <- arch::arch_allocate_schema()
    array <- arch::arch_allocate_array_data()
    ## modeled after libtiledb_query_export_buffer_arch_pointers
    export_column(uri, column, schema, array)
    aa <- arch::arch_array(schema, array, FALSE)
    res <- arch::from_arch_array(aa)
    res
}

##' @rdname get_table
arrow_to_dt <- function(arr) {
    data.table::data.table(as.data.frame(dplyr::collect(arrow::as_arrow_table(arch::from_arch_array(arr, arrow::RecordBatch)))))
}

##' @importFrom Rcpp evalCpp
##' @importFrom arch arch_allocate_schema arch_allocate_array_data arch_array as_arch_array_stream from_arch_array arch_schema_info
##' @importFrom data.table data.table
##' @importFrom dplyr collect
##' @exportPattern "^[[:alpha:]]+"
##' @useDynLib tiledbsoma, .registration=TRUE
NULL
