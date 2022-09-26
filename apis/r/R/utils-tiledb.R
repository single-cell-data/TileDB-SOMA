tiledb_zstd_filter <- function(level = 3L) {
  tiledb::tiledb_filter_set_option(
    object = tiledb::tiledb_filter("ZSTD"),
    option = "COMPRESSION_LEVEL",
    value = level
  )
}
