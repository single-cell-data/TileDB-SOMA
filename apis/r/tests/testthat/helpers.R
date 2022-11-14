# Convert all factor columns in a data.frame to characters
fac2char <- function(x) {
  stopifnot(is.data.frame(x))
  factcols <- vapply_lgl(x, is.factor)
  x[factcols] <- lapply(x[factcols], as.character)
  return(x)
}


# Temporarily set the allocation size preference
with_allocation_size_preference <- function(value, .local_envir = parent.frame()) {
  orig_value <- tiledb::get_allocation_size_preference()
  withr::defer(
    tiledb::set_allocation_size_preference(orig_value),
    envir = .local_envir
  )
  tiledb::set_allocation_size_preference(value)
}


#' Compare data.frames or matrices after sorting dimensions
#' Useful for comparing objects were numerically or lexicographically sorted
#' after being reconstituted from TileDB.
#' @noRd
expect_equal_after_ordering_dims <- function(object, expected, ...) {
  stopifnot(
    is.data.frame(object) || is_matrix(object),
    is.data.frame(expected) || is_matrix(expected)
  )

  object <- object[rownames(expected), colnames(expected)]
  testthat::expect_equal(object, expected, ...)
}
