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

check_package <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible())
  }
  stop(paste0("Package '", package, "' must be installed"))
}

#' Assert all values of `x` are a subset of `y`. @param x,y vectors of values
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
