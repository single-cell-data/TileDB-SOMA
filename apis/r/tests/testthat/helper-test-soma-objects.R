create_and_populate_soma_dataframe <- function(
  uri,
  nrows = 10L,
  seed = 1,
  index_column_names = "foo"
) {
  set.seed(seed)

  arrow_schema <- create_arrow_schema()

  tbl <- arrow::arrow_table(
    foo = seq.int(nrows) + 1000L,
    soma_joinid = bit64::seq.integer64(from = 0L, to = nrows - 1L),
    bar = seq(nrows) + 0.1,
    baz = as.character(seq.int(nrows) + 1000L),
    schema = arrow_schema
  )

  sdf <- SOMADataFrame$new(uri, internal_use_only = "allowed_use")
  sdf$create(arrow_schema, index_column_names = index_column_names)
  sdf$write(tbl)
  sdf
}

create_and_populate_obs <- function(uri, nrows = 10L, seed = 1) {
  create_and_populate_soma_dataframe(
    uri = uri,
    nrows = nrows,
    seed = seed,
    index_column_names = "soma_joinid"
  )
}

create_and_populate_var <- function(uri, nrows = 10L, seed = 1) {

  tbl <- arrow::arrow_table(
    soma_joinid = bit64::seq.integer64(from = 0L, to = nrows - 1L),
    quux = as.character(seq.int(nrows) + 1000L),
    xyzzy = runif(nrows),
    schema = arrow::schema(
      arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
      arrow::field("quux", arrow::large_utf8(), nullable = FALSE),
      arrow::field("xyzzy", arrow::float64(), nullable = FALSE)
    )
  )

  sdf <- SOMADataFrame$new(uri, internal_use_only = "allowed_use")
  sdf$create(tbl$schema, index_column_names = "soma_joinid")
  sdf$write(tbl)
  sdf
}

#' @param ... Arguments passed to create_sparse_matrix_with_int_dims
create_and_populate_sparse_nd_array <- function(uri, ...) {

  smat <- create_sparse_matrix_with_int_dims(...)

  ndarray <- SOMASparseNDArray$new(uri, internal_use_only = "allowed_use")
  ndarray$create(arrow::int32(), shape = dim(smat))
  ndarray$write(smat)
  ndarray
}

# Create a SOMAExperiment with a single measurement, "RNA"
# Example with X_layer_names = c("counts", "logcounts"):
#  soma-experiment-query-all1c20a1d341584 GROUP
#  |-- obs ARRAY
#  |-- ms GROUP
#  |------ RNA GROUP
#  |---------- var ARRAY
#  |---------- X GROUP
#  |-------------- counts ARRAY
#  |-------------- logcounts ARRAY
create_and_populate_experiment <- function(
  uri,
  n_obs,
  n_var,
  X_layer_names,
  config = NULL
) {

  experiment <- SOMAExperiment$new(uri, platform_config = config, internal_use_only = "allowed_use")$create()
  experiment$obs <- create_and_populate_obs(
    uri = file.path(uri, "obs"),
    nrows = n_obs
  )
  experiment$ms <- SOMACollection$new(file.path(uri, "ms"), internal_use_only = "allowed_use")$create()

  ms_rna <- SOMAMeasurement$new(file.path(uri, "ms", "RNA"), internal_use_only = "allowed_use")$create()
  ms_rna$var <- create_and_populate_var(
    uri = file.path(ms_rna$uri, "var"),
    nrows = n_var
  )
  ms_rna$X <- SOMACollection$new(file.path(ms_rna$uri, "X"), internal_use_only = "allowed_use")$create()

  for (layer_name in X_layer_names) {
    nda <- create_and_populate_sparse_nd_array(
      uri = file.path(ms_rna$X$uri, layer_name),
      nrows = n_obs,
      ncols = n_var
    )
    ms_rna$X$set(nda, name = layer_name)
  }

  experiment$ms$set(ms_rna, name = "RNA")
  experiment
}
