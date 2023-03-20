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

#' Pad a sparse Matrix with additional rows or columns
#'
#' @param x A dgTMatrix
#' @param colnames,rownames A vector of column or row names
#' to add to the matrix.
#' @param returns A padded matrix containing all provided
#' row/column names
#'
#' @importFrom Matrix sparseMatrix
#'
#' @noRd
#'
pad_matrix <- function(x, rownames = NULL, colnames = NULL) {
  stopifnot(
    inherits(x, "Matrix"),
    is.character(colnames) || is.character(rownames)
  )
  # lookup table for Matrix representations
  mat_rep <- switch(
    EXPR = class(x),
    dgTMatrix = "T",
    dgCMatrix = "C",
    dgRMatrix = "R",
    stop("Untested Matrix object representation")

  )
  new_rownames <- setdiff(rownames, rownames(x))
  new_colnames <- setdiff(colnames, colnames(x))
  dtype <- typeof(methods::slot(object = x, name = 'x'))
  if (!is_empty(new_rownames)) {
    rpad <- Matrix::sparseMatrix(
      i = integer(0L),
      j = integer(0L),
      x = vector(mode = dtype, length = 0L),
      dims = c(length(new_rownames), ncol(x)),
      dimnames = list(new_rownames, colnames(x)),
      repr = mat_rep
    )
    x <- rbind(x, rpad)
  }
  if (!is_empty(new_colnames)) {
    cpad <- Matrix::sparseMatrix(
      i = integer(0L),
      j = integer(0L),
      x = vector(mode = dtype, length = 0L),
      dims = c(nrow(x), length(new_colnames)),
      dimnames = list(rownames(x), new_colnames),
      repr = mat_rep
    )
    x <- cbind(x, cpad)
  }
  x
}

#' Get R Version
#'
#' @param repr Representation of R version; choose from:
#' \itemize{
#'  \item \dQuote{\code{v}}: a \code{\link[base]{package_version}}
#'  \item \dQuote{\code{c}}: a character
#' }
#'
#' @return The version of R currently being used
#'
#' @keywords internal
#'
#' @noRd
#'
r_version <- function(repr = c('v', 'c')) {
  repr <- repr[1L]
  repr <- match.arg(arg = repr)
  version <- paste(R.Version()[c('major', 'minor')], collapse = '.')
  return(switch(EXPR = repr, v = package_version(version), version))
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
