create_and_populate_obs <- function(uri) {
  create_and_populate_soma_dataframe(uri)
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
