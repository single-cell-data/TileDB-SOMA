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

err_to_warn <- function(err, immediate. = TRUE) {
  warning(conditionMessage(err), call. = FALSE, immediate. = immediate.)
  return(invisible(err))
}

null <- function(...) {
  return(NULL)
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
SOMA_ENCODING_VERSION <- "1"

check_arrow_pointers <- function(arrlst) {
    stopifnot("First argument must be an external pointer to ArrowArray" = check_arrow_array_tag(arrlst[[1]]),
              "Second argument must be an external pointer to ArrowSchema" = check_arrow_schema_tag(arrlst[[2]]))
}

#' @importFrom Matrix as.matrix
#' @importFrom arrow RecordBatch
#' @import R6 methods utils
##' @importFrom Rcpp evalCpp
##' @importFrom data.table data.table
##' @importFrom spdl setup
##' @useDynLib tiledbsoma, .registration=TRUE
NULL
