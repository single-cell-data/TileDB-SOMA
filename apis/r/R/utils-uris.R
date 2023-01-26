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
