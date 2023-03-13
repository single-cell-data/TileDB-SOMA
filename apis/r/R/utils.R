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
is_named <- function(x) {
  !is.null(names(x))
}

is_named_list <- function(x) {
  is.list(x) && is_named(x)
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

# a matrix or one of the Matrix package matrices
is_matrix <- function(x) {
  is.matrix(x) || inherits(x, "Matrix")
}

# a matrix/Matrix with non-empty dimension names
is_labeled_matrix <- function(x) {
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

check_package <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible())
  }
  stop(paste0("Package '", package, "' must be installed"))
}

is_remote_uri <- function(x) {
  string_starts_with(x, "s3://") | string_starts_with(x, "tiledb://")
}

# Drop-in replacement for file.paths() that ignores the platform separator when
# constructing remote S3 or TileDB URIs
file_path <- function(..., fsep = .Platform$file.sep) {
  paths <- list(...)
  if (is_remote_uri(paths[[1]])) fsep <- "/"
  file.path(..., fsep = fsep)
}

#' Assert all values of `x` are a subset of `y`.
#' @param x,y vectors of values
#' @param type A character vector of length 1 used in the error message
#' @return `TRUE` if all values of `x` are present in `y`, otherwise an
#' informative error is thrown with the missing values.
#' @noRd
assert_subset <- function(x, y, type = "value") {
  stopifnot(is.atomic(x) && is.atomic(y))
  missing <- !x %in% y
  if (any(missing)) {
    stop(sprintf(
      "The following %s%s not exist: %s",
      type,
      ifelse(length(missing) == 1, " does", "s do"),
      glue::glue_collapse(x[missing], sep = ", ", last = " and ")
    ), call. = FALSE)
  }
  TRUE
}

SOMA_OBJECT_TYPE_METADATA_KEY <- "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY <- "soma_encoding_version"
SOMA_ENCODING_VERSION <- "0"
SOMA_LEGACY_VALIDITY_KEY <- "soma_legacy_validity"
SOMA_LEGACY_VALIDITY <- "false"
TILEDB_LEGACY_KEY <- "r.legacy_validity_mode"
