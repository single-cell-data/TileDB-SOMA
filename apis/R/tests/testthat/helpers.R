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
