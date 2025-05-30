tiledb_zstd_filter <- function(level = 3L) {
  tiledb::tiledb_filter_set_option(
    object = tiledb::tiledb_filter("ZSTD"),
    option = "COMPRESSION_LEVEL",
    value = level
  )
}

# Notes:
# * "auto" is used by libtiledbsoma -- it's not a core-level tiledb parameter
# * core-level "GLOBAL_ORDER" and "UNORDERED" are not currently supported by libtiledbsoma
match_query_layout <- function(layout) {
  layouts <- c("ROW_MAJOR", "COL_MAJOR", "auto")
  match.arg(layout, choices = layouts, several.ok = FALSE)
}

map_query_layout <- function(layout) {
  switch(layout,
    ROW_MAJOR = "row-major",
    COL_MAJOR = "column-major",
    tolower(layout)
  )
}

get_tiledb_version <- function(compact = FALSE) {
  stopifnot("'compact' must be TRUE or FALSE" = isTRUE(compact) || isFALSE(compact))
  version <- `names<-`(
    tiledb_embedded_version(),
    c("major", "minor", "patch")
  )
  if (compact) {
    return(paste(version, collapse = "."))
  }
  return(version)
}

#' Display package versions
#'
#' Print version information for \CRANpkg{tiledb} (R package), libtiledbsoma,
#' and TileDB embedded, suitable for assisting with bug reports.
#'
#' @export
#'
#' @examples
#' show_package_versions()
#'
show_package_versions <- function() {
  cat("tiledbsoma:    ", toString(utils::packageVersion("tiledbsoma")), "\n",
    "tiledb-r:      ", toString(utils::packageVersion("tiledb")), "\n",
    "tiledb core:   ", as.character(get_tiledb_version(compact = TRUE)), "\n",
    "libtiledbsoma: ", libtiledbsoma_version(compact = TRUE), "\n",
    "R:             ", R.version.string, "\n",
    "OS:            ", utils::osVersion, "\n",
    sep = ""
  )
}

#' @rdname tiledbsoma_stats
#'
#' @export
#'
tiledbsoma_stats_show <- function() {
  cat(tiledbsoma_stats_dump(), "\n")
}
