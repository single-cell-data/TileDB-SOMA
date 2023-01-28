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

# Create a SOMAExperiment with a single measurement, "RNA"
create_and_populate_experiment <- function(uri, n_obs, n_var, X_layer_names) {

  experiment <- SOMAExperiment$new(uri)$create()
  experiment$obs <- create_and_populate_obs(
    uri = file.path(uri, "obs"),
    nrows = n_obs
  )
  experiment$ms <- SOMACollection$new(file.path(uri, "ms"))$create()

  ms_rna <- SOMAMeasurement$new(file.path(uri, "ms", "RNA"))$create()
  ms_rna$var <- create_and_populate_var(file.path(uri, "var"), nrows = n_var)
  ms_rna$X <- SOMACollection$new(file.path(uri, "ms", "X"))$create()

  for (layer_name in X_layer_names) {
    nda <- create_and_populate_sparse_nd_array(
      uri = file.path(uri, "ms", "X", layer_name),
      nrows = n_obs,
      ncols = n_var
    )
    ms_rna$X$set(nda, name = layer_name)
  }

  experiment$ms$set(ms_rna, name = "RNA")
  experiment
}
