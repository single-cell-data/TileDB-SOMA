##' @noRd
get_table <- function(uri) {
    colnames <- get_column_names(uri)
    ll <- lapply(colnames, function(n) get_column(uri, n))
    names(ll) <- colnames
    as.data.frame(ll)
}

##' @noRd
get_column <- function(uri, column) {
    schema <- arch::arch_allocate_schema()
    array <- arch::arch_allocate_array_data()
    ## modeled after libtiledb_query_export_buffer_arch_pointers
    export_column(uri, column, schema, array)
    aa <- arch::arch_array(schema, array, FALSE)
    res <- arch::from_arch_array(aa)
    res
}

##' @importFrom Rcpp evalCpp
##' @importFrom arch arch_allocate_schema arch_allocate_array_data arch_array as_arch_array_stream from_arch_array arch_schema_info
##' @useDynLib tiledbsoma, .registration=TRUE
NULL
