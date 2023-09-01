create_empty_test_array <- function(uri) {
  stopifnot(!dir.exists(uri))
  dim <- tiledb::tiledb_dim("d0", type = "ASCII", domain = NULL, tile = NULL)
  dom <- tiledb::tiledb_domain(dims = dim)
  schema <- tiledb::tiledb_array_schema(
    domain = dom,
    attrs = c(tiledb::tiledb_attr("a", type = "INT32")),
    sparse = TRUE
  )
  tiledb::tiledb_array_create(uri, schema)
  return(uri)
}

extended_tests <- function() {
    ## check if at CI, if so extended test
    ## could add if pre-release number ie 1.4.3.1 instead of 1.4.3
    Sys.getenv("CI", "") != ""
}
