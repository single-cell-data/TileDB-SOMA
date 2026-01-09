is_remote_uri <- function(x) {
  string_starts_with(x, "s3://") | string_starts_with(x, "tiledb://")
}

# Drop-in replacement for file.paths() that ignores the platform separator when
# constructing remote S3 or TileDB URIs
file_path <- function(..., fsep = .Platform$file.sep) {
  paths <- list(...)
  if (is_remote_uri(paths[[1]])) {
    fsep <- "/"
  }
  file.path(..., fsep = fsep)
}

#' Return the scheme of a URI
#' @noRd
uri_scheme <- function(uri) {
  stopifnot(is_scalar_character(uri))
  uri_parts <- strsplit(uri, "://")[[1]]
  if (length(uri_parts) == 1) {
    return(NULL)
  }
  uri_parts[[1]]
}

#' Remove the scheme from a URI
#' @noRd
uri_scheme_remove <- function(uri) {
  stopifnot(is_scalar_character(uri))
  uri_parts <- strsplit(uri, "://")[[1]]
  if (length(uri_parts) == 1) {
    return(uri)
  }
  uri_parts[[2]]
}

#' Return a URI relative to a parent URI
#' This takes URI schemes into account and errors if they do not match. URIs
#' without a scheme are treated as `file://` URIs.
#'
#' @noRd
#'
make_uri_relative <- function(uri, relative_to) {
  stopifnot(
    "'uri' and 'relative_to' must be scalar character vectors" = is_scalar_character(
      uri
    ) &&
      is_scalar_character(relative_to)
  )

  uri_scheme <- uri_scheme(uri)
  relative_to_scheme <- uri_scheme(relative_to)
  if (uri_scheme %||% "file" != relative_to_scheme %||% "file") {
    stop(
      "Unable to make relative path between URIs with different schemes",
      call. = FALSE
    )
  }

  # Remove schemes from URIs before calculating relative path
  uri <- uri_scheme_remove(uri)
  relative_to <- uri_scheme_remove(relative_to)

  if (!fs::path_has_parent(uri, relative_to)) {
    stop(
      "Unable to make relative path between URIs with no common parent",
      call. = FALSE
    )
  }

  fs::path_rel(
    path = uri_scheme_remove(uri),
    start = uri_scheme_remove(relative_to)
  )
}

#' Detects whether a provided URI is a child relative URI to the parent.
#' @noRd
is_relative_uri <- function(uri) {
  stopifnot(is_scalar_character(uri))
  return(
    !grepl("://", uri, fixed = TRUE) &&
      !startsWith(uri, "/") &&
      !startsWith(uri, "~")
  )
}
