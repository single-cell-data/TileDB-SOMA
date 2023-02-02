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

#' show_package_versions
#'
#' @description Prints package information suitable for assisting with bug reports.
#'
#' @export
show_package_versions <- function() {
    cat("tiledbsoma:   ", toString(utils::packageVersion("tiledbsoma")),    "\n")
    cat("tiledb:       ", toString(utils::packageVersion("tiledb")),      "\n")
    cat("R:            ", R.Version()$version.string,                     "\n")
}
