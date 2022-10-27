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

##' @noRd
arrow_to_dt <- function(arr) {
    data.table::data.table(as.data.frame(dplyr::collect(arrow::as_arrow_table(arch::from_arch_array(arr, arrow::RecordBatch)))))
}

##' @importFrom Rcpp evalCpp
##' @importFrom arch arch_allocate_schema arch_allocate_array_data arch_array as_arch_array_stream from_arch_array arch_schema_info
##' @importFrom data.table data.table
##' @importFrom dplyr collect
##' @useDynLib tiledbsoma, .registration=TRUE
NULL
