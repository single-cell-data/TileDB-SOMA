##' Access a DataFrame or Single Column From a Given URI
##'
##' These functions access a given SOMA URI and extract, respectively, the
##' complete data.frame or named column from the data set.
##'
##' @param uri Character value with URI path to a SOMA data set
##' @param column Character value with the name of the column to retrieve
##' @return The selected column from the give data set
##' @examples
##' \dontrun{
##' uri <- "test/soco/pbmc3k_processed/obs"
##' column <-  "n_counts"
##' summary(getColumn(uri, column)
##' summary(getTable(uri))
##' }
##' @export
getColumn <- function(uri, column) {
    schema <- arch::arch_allocate_schema()
    array <- arch::arch_allocate_array_data()
    ## modeled after libtiledb_query_export_buffer_arch_pointers
    export_column(uri, column, schema, array)
    aa <- arch::arch_array(schema, array, FALSE)
    res <- arch::from_arch_array(aa)
    res
}

##' @rdname getColumn
##' @export
getTable <- function(uri = "test/soco/pbmc3k_processed/obs") {
    colnames <- get_column_names(uri)
    ll <- lapply(colnames, \(n) getColumn(uri, n))
    names(ll) <- colnames
    as.data.frame(ll)
}

##' @importFrom Rcpp evalCpp
##' @useDynLib tiledbsoma, .registration=TRUE
NULL
