#' Check if object is empty
#' @noRd
is_empty <- function(x) {
  switch(class(x)[1],
    "data.frame" = nrow(x) == 0,
    length(x) == 0
  )
}

#' Check if a vector is named
#' @noRd
is_named <- function(x, allow_empty = TRUE) {
  !is.null(names(x)) && ifelse(allow_empty, TRUE, all(nzchar(x = names(x = x))))
}

is_named_list <- function(x) {
  is.list(x) && is_named(x)
}

is_scalar_logical <- function(x) {
  is.logical(x) && length(x) == 1
}

is_scalar_character <- function(x) {
  is.character(x) && length(x) == 1
}

is_character_or_null <- function(x) {
  is.character(x) || is.null(x)
}

has_character_rownames <- function(x) {
  stopifnot(is.data.frame(x))
  typeof(attr(x, "row.names")) == "character"
}

is_matrix <- function(x) {
  is.matrix(x) || inherits(x, "Matrix")
}

is_vector_or_int64 <- function(x) {
    is.vector(x) || inherits(x, "integer64")
}

has_dimnames <- function(x) {
  stopifnot(is_matrix(x))
  dims <- dimnames(x) %||% list(NULL)
  all(!vapply(dims, is.null, logical(1L)))
}

string_starts_with <- function(x, prefix) {
  prefix <- paste0("^", prefix)
  grepl(prefix, x)
}

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

check_arrow_pointers <- function(arrlst) {
    stopifnot("First argument must be an external pointer to ArrowArray" = check_arrow_array_tag(arrlst[[1]]),
              "Second argument must be an external pointer to ArrowSchema" = check_arrow_schema_tag(arrlst[[2]]))
}

##' @noRd
arrow_to_dt <- function(arrlst) {
    check_arrow_pointers(arrlst)
    rb <- dplyr::collect(arrow::RecordBatch$import_from_c(arrlst[[1]], arrlst[[2]]))
    data.table(as.data.frame(rb))
}

##' @noRd
as_arrow_table <- function(arrlst) {
    check_arrow_pointers(arrlst)
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
