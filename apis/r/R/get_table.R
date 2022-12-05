##' @rdname soma_reader
arrow_to_dt <- function(arrlst) {
    ## this helper will be replaced once the under-development package 'nanoarrow' (on
    ## github at apache/arror-nanoarrow) is released, for now we use 'arch' which predates it
    data.table::data.table(dplyr::collect(arch::from_arch_array(arrlst, arrow::RecordBatch)))
}

##' @importFrom Rcpp evalCpp
##' @importFrom arch arch_allocate_schema arch_allocate_array_data arch_array as_arch_array_stream from_arch_array arch_schema_info
##' @importFrom data.table data.table
##' @importFrom dplyr collect
##' @importFrom spdl setup
##' @exportPattern "^[[:alpha:]]+"
##' @useDynLib tiledbsoma, .registration=TRUE
NULL
