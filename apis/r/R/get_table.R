##' Access a DataFrame or Single Column From a Given URI
##'
##' These functions access a given SOMA URI and extract, respectively, the
##' complete data.frame or named column from the data set.
##'
##' @param uri Character value with URI path to a SOMA data set
##' @param column Character value with the name of the column to retrieve
##' @param colnames Vector of character value with the name of the columns to retrieve
##' @param qc External Pointer object to TileDB Query Condition, defaults to \sQuote{NULL} i.e.
##' no query condition
##' @param loglevel Character value with the desired logging level, defaults to \sQuote{warn}
##' @return The selected data object(s) from the given data set
##' @examples
##' \dontrun{
##' uri <- "test/soco/pbmc3k_processed/obs"
##' column <-  "n_counts"
##' summary(get_table(uri))
##' summary(get_columns(uri, column))
##' columns <- c("n_genes", "louvain")
##' z <- export_recordbatch(uri, columns)
##' rb <- arch::from_arch_array(z, arrow::RecordBatch)
##' z <- export_recordbatch(uri, columns)
##' tb <- arrow::as_arrow_table(arch::from_arch_array(z, arrow::RecordBatch))
##' }
##' @importFrom arch arch_allocate_schema arch_allocate_array_data arch_array as_arch_array_stream from_arch_array arch_schema_info
##' @export
get_table <- function(uri) {
    colnames <- get_column_names(uri)
    ll <- lapply(colnames, function(n) get_column(uri, n))
    names(ll) <- colnames
    as.data.frame(ll)
}

##' @rdname get_table
##' @export
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
##' @useDynLib tiledbsoma, .registration=TRUE
NULL
