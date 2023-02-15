tiledb_zstd_filter <- function(level = 3L) {
  tiledb::tiledb_filter_set_option(
    object = tiledb::tiledb_filter("ZSTD"),
    option = "COMPRESSION_LEVEL",
    value = level
  )
}

match_query_layout <- function(layout) {
  layouts <- c("ROW_MAJOR", "COL_MAJOR", "GLOBAL_ORDER", "UNORDERED")
  match.arg(layout, choices = layouts, several.ok = FALSE)
}

map_query_layout <- function(layout) {
    switch(layout,
           ROW_MAJOR = "row-major",
           COL_MAJOR = "column-major",
           tolower(layout)
    )
}

#' Display Package Versions
#'
#' This helperfunction prints package information suitable for assisting with bug reports.
#'
#' @export
#' @importFrom utils packageVersion
show_package_versions <- function() {
    cat("tiledbsoma:   ", toString(utils::packageVersion("tiledbsoma")), "\n")
    cat("tiledb-r:     ", toString(utils::packageVersion("tiledb")), "\n")
    cat("tiledb core:  ", as.character(tiledb::tiledb_version(compact=TRUE)), "\n")
    cat("R:            ", R.version.string, "\n")
}

#' @rdname tiledbsoma_stats_enable
#' @export
tiledbsoma_stats_show <- function() {
    cat(stats_dump(), "\n")
}
