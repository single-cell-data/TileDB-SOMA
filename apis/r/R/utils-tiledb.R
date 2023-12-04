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

#' Display package versions
#'
#' Print version information for \pkg{tiledb} (R package), libtiledbsoma, and
#' TileDB embedded, suitable for assisting with bug reports.
#'
#' @export
#' @importFrom utils packageVersion
show_package_versions <- function() {
    cat("tiledbsoma:    ", toString(utils::packageVersion("tiledbsoma")), "\n",
        "tiledb-r:      ", toString(utils::packageVersion("tiledb")), "\n",
        "tiledb core:   ", as.character(tiledb::tiledb_version(compact=TRUE)), "\n",
        "libtiledbsoma: ", libtiledbsoma_version(compact=TRUE), "\n",
        "R:             ", R.version.string, "\n",
        "OS:            ", utils::osVersion, "\n",
        sep="")
}

#' @rdname tiledbsoma_stats
#' @export
tiledbsoma_stats_show <- function() {
    cat(tiledbsoma_stats_dump(), "\n")
}
