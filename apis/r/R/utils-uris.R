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

#' Return the scheme of a URI
#' @noRd
uri_scheme <- function(uri) {
  stopifnot(is_scalar_character(uri))
  uri_parts <- strsplit(uri, "://")[[1]]
  if (length(uri_parts) == 1) return(NULL)
  uri_parts[[1]]
}

#' Remove the scheme from a URI
#' @noRd
uri_scheme_remove <- function(uri) {
  stopifnot(is_scalar_character(uri))
  uri_parts <- strsplit(uri, "://")[[1]]
  if (length(uri_parts) == 1) return(uri)
  uri_parts[[2]]
}
