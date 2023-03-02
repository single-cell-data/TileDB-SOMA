#' @importFrom glue glue_collapse
string_collapse <- function(x, sep = ", ") {
  glue::glue_collapse(x, sep = ", ", width = getOption("width", Inf))
}

n_unique <- function(x) {
  length(unique(x))
}

vapply_char <- function(X, FUN, ..., USE.NAMES = TRUE) {
  vapply(X, FUN, FUN.VALUE = character(1L), ..., USE.NAMES = USE.NAMES)
}

vapply_lgl <- function(X, FUN, ..., USE.NAMES = TRUE) {
  vapply(X, FUN, FUN.VALUE = logical(1L), ..., USE.NAMES = USE.NAMES)
}

vapply_int <- function(X, FUN, ..., USE.NAMES = TRUE) {
  vapply(X, FUN, FUN.VALUE = integer(1L), ..., USE.NAMES = USE.NAMES)
}

# rename(iris, c(petal_length = "Petal.Length", species = "Species", hi = "YO"))
rename <- function(x, names) {
  stopifnot(
    "'x' must be named" = is_named(x),
    "'names' must be a named character vector" = is_named(names),
    "All 'names' must be in 'x'" = all(names %in% names(x))
  )

  name_index <- match(names, names(x))
  names(x)[name_index] <- names(names)
  x
}

# Return y if x is NULL, else x
`%||%` <- function(x, y) {
  if (missing(x) || is.null(x) || length(x) == 0) y else x
}

# For use in read-only R6 active bindings
read_only_error <- function(field_name) {
  stop(
    sprintf("'%s' is a read-only field.", field_name),
    call. = FALSE
  )
}

SOMA_OBJECT_TYPE_METADATA_KEY <- "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY <- "soma_encoding_version"
SOMA_ENCODING_VERSION <- "0"

##' @rdname soma_reader
arrow_to_dt <- function(arrlst) {
    ## this helper will be replaced once the under-development package 'nanoarrow' (on
    ## github at apache/arror-nanoarrow) is released, for now we use 'arch' which predates it
    #data.table::data.table(dplyr::collect(arch::from_arch_array(arrlst, arrow::RecordBatch)))
    #data.table::data.table(dplyr::collect(arch::from_arch_array(arrlst, arrow::RecordBatch)))
    rb <- dplyr::collect(arrow::RecordBatch$import_from_c(arrlst[[1]], arrlst[[2]]))
    #data.table::data.table(dplyr::collect()) #arch::from_arch_array(arrlst, arrow::RecordBatch)))
    data.table(as.data.frame(rb))
}

##' @rdname soma_reader
as_arrow_table <- function(arrlst) {
    arrow::as_arrow_table(arrow::RecordBatch$import_from_c(arrlst[[1]], arrlst[[2]]))
}

#' @importFrom Matrix as.matrix
#' @importFrom arrow RecordBatch
#' @import R6 methods utils
##' @importFrom Rcpp evalCpp
##' @importFrom data.table data.table
##' @importFrom dplyr collect
##' @importFrom spdl setup
##' @useDynLib tiledbsoma, .registration=TRUE
NULL
