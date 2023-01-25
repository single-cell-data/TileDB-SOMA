create_and_populate_obs <- function(uri, nrows = 10L, seed = 1) {
  create_and_populate_soma_dataframe(uri, nrows, seed)
}

create_and_populate_var <- function(uri, nrows = 10L, seed = 1) {

  tbl <- arrow::arrow_table(
    soma_joinid = bit64::seq.integer64(from = 1L, to = nrows),
    quux = sample(letters, size = nrows, replace = TRUE),
    xyzzy = runif(nrows)
  )

  sdf <- SOMADataFrame$new(uri)
  sdf$create(tbl$schema, index_column_names = "soma_joinid")
  sdf$write(tbl)
  sdf
}

create_and_populate_var <- function(uri) {

  tbl <- arrow::arrow_table(
    soma_joinid = bit64::as.integer64(1L:4L),
    quux = c("alvin", "simon", "theodore", "dave"),
    xyzzy = c(12.3, 45.6, 78.9, 10.11)
  )

  sdf <- SOMADataFrame$new(uri)
  sdf$create(tbl$schema, index_column_names = "soma_joinid")
  sdf$write(tbl)
  sdf
}
