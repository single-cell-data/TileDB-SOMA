# Returns the object created, populated, and closed (unless otherwise requested)
create_and_populate_soma_dataframe <- function(
  uri,
  nrows = 10L,
  seed = 1,
  index_column_names = "foo",
  mode = NULL
) {
  set.seed(seed)

  arrow_schema <- create_arrow_schema()
  tbl <- create_arrow_table(nrows = nrows)

  sdf <- SOMADataFrameCreate(uri, arrow_schema, index_column_names = index_column_names)
  sdf$write(tbl)

  if (is.null(mode)) {
    sdf$close()
  } else if (mode == "READ") {
    sdf$close()
    sdf <- SOMADataFrameOpen(uri, mode = mode)
  }
  sdf
}

# Returns the object created, populated, and closed (unless otherwise requested)
create_and_populate_obs <- function(uri, nrows = 10L, seed = 1, mode = NULL) {
  create_and_populate_soma_dataframe(
    uri = uri,
    nrows = nrows,
    seed = seed,
    index_column_names = "soma_joinid"
  )
}

# Returns the object created, populated, and closed (unless otherwise requested)
create_and_populate_var <- function(uri, nrows = 10L, seed = 1, mode = NULL) {

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

  sdf <- SOMADataFrameCreate(uri, tbl$schema, index_column_names = "soma_joinid")
  sdf$write(tbl)

  if (is.null(mode)) {
    sdf$close()
  } else if (mode == "READ") {
    sdf$close()
    sdf <- SOMADataFrameOpen(uri, mode = mode)
  }
  sdf
}

# Creates a SOMAExperiment with a single measurement, "RNA"
# Returns the object created, populated, and closed (unless otherwise requested)
#' @param ... Arguments passed to create_sparse_matrix_with_int_dims
create_and_populate_sparse_nd_array <- function(uri, mode = NULL, ...) {
  smat <- create_sparse_matrix_with_int_dims(...)

  ndarray <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = dim(smat))
  ndarray$write(smat)

  if (is.null(mode)) {
    ndarray$close()
  } else if (mode == "READ") {
    ndarray$close()
    ndarray <- SOMASparseNDArrayOpen(uri, mode = mode)
  }
  ndarray
}

# Creates a SOMAExperiment with a single measurement, "RNA"; populates it;
# returns it closed (unless otherwise requested).
#
# Example with X_layer_names = c("counts", "logcounts"):
#  soma-experiment-query-all1c20a1d341584 GROUP
#  |-- obs ARRAY
#  |-- ms GROUP
#  |------ RNA GROUP
#  |---------- var ARRAY
#  |---------- X GROUP
#  |-------------- counts ARRAY
#  |-------------- logcounts ARRAY
#' @param obsm_layers A named integer vector of layers to add to `obsm`; the
#' names will be used to create new layers in `obsm` and the value will determine
#' the number of dimensions (columns) to add. Names starting with `dense:` will
#' be created as _dense_ arrays (eg. `dense:X_ica`). Pass `NULL` to prevent
#' creation of `obsm` layers
#' @param varm_layers A named integer vector of layers to add to `varm`; the
#' names will be used to create new layers in `varm` and the value will determine
#' the number of dimensions (columns) to add. Names starting with `dense:` will
#' be created as _dense_ arrays (eg. `dense:ICs`). Pass `NULL` to prevent
#' creation of `varm` layers
#' @param obsp_layers A character vector of `obsp` layers; pass `NULL` to
#' prevent creation of `obsp` layers
#' @param varp_layers A character vector of `varp` layers; pass `NULL` to
#' prevent creation of `varp` layers
#'
create_and_populate_experiment <- function(
  uri,
  n_obs,
  n_var,
  X_layer_names,
  obsm_layers = NULL,
  varm_layers = NULL,
  obsp_layer_names = NULL,
  varp_layer_names = NULL,
  config = NULL,
  mode = NULL
) {

  stopifnot(
    "'obsm_layers' must be a named integer vector" = is.null(obsm_layers) ||
      (rlang::is_integerish(obsm_layers) && rlang::is_named(obsm_layers) && all(obsm_layers > 0L)),
    "'varm_layers' must be a named integer vector" = is.null(varm_layers) ||
      (rlang::is_integerish(varm_layers) && rlang::is_named(varm_layers) && all(varm_layers > 0L)),
    "'obsp_layer_names' must be a character vector" = is.null(obsp_layer_names) ||
      (is.character(obsp_layer_names) && all(nzchar(obsp_layer_names))),
    "'varp_layer_names' must be a character vector" = is.null(varp_layer_names) ||
      (is.character(varp_layer_names) && all(nzchar(varp_layer_names)))
  )

  experiment <- SOMAExperimentCreate(uri, platform_config = config)

  experiment$obs <- create_and_populate_obs(
    uri = file.path(uri, "obs"),
    nrows = n_obs
  )

  experiment$ms <- SOMACollectionCreate(file.path(uri, "ms"))

  ms_rna <- SOMAMeasurementCreate(file.path(uri, "ms", "RNA"))
  ms_rna$var <- create_and_populate_var(
    uri = file.path(ms_rna$uri, "var"),
    nrows = n_var
  )
  ms_rna$X <- SOMACollectionCreate(file.path(ms_rna$uri, "X"))

  for (layer_name in X_layer_names) {
    snda <- create_and_populate_sparse_nd_array(
      uri = file.path(ms_rna$X$uri, layer_name),
      nrows = n_obs,
      ncols = n_var
    )
    ms_rna$X$set(snda, name = layer_name)
  }
  ms_rna$X$close()

  # Add obsm layers
  if (rlang::is_integerish(obsm_layers)) {
    obsm <- SOMACollectionCreate(file.path(ms_rna$uri, "obsm"))
    for (layer in names(obsm_layers)) {
      key <- gsub(pattern = '^dense:', replacement = '', x = layer)
      shape <- c(n_obs, obsm_layers[layer])
      if (grepl(pattern = '^dense:', x = layer)) {
        obsm$add_new_dense_ndarray(
          key = key,
          type = arrow::int32(),
          shape = shape
        )
        obsm$get(key)$write(create_dense_matrix_with_int_dims(
          nrows = shape[1L],
          ncols = shape[2L]
        ))
      } else {
        obsm$add_new_sparse_ndarray(
          key = key,
          type = arrow::int32(),
          shape = shape
        )
        obsm$get(key)$write(create_sparse_matrix_with_int_dims(
          nrows = shape[1L],
          ncols = shape[2L]
        ))
      }
    }
    obsm$close()
    ms_rna$add_new_collection(obsm, "obsm")
  }

  # Add varm layers
  if (rlang::is_integerish(varm_layers)) {
    varm <- SOMACollectionCreate(file.path(ms_rna$uri, "varm"))
    for (layer in names(varm_layers)) {
      key <- gsub(pattern = '^dense:', replacement = '', x = layer)
      shape <- c(n_var, varm_layers[layer])
      if (grepl(pattern = '^dense:', x = layer)) {
        varm$add_new_dense_ndarray(
          key = key,
          type = arrow::int32(),
          shape = shape
        )
        varm$get(key)$write(create_dense_matrix_with_int_dims(
          nrows = shape[1L],
          ncols = shape[2L]
        ))
      } else {
        varm$add_new_sparse_ndarray(
          key = key,
          type = arrow::int32(),
          shape = shape
        )
        varm$get(key)$write(create_sparse_matrix_with_int_dims(
          nrows = shape[1L],
          ncols = shape[2L]
        ))
      }
    }
    varm$close()
    ms_rna$add_new_collection(varm, "varm")
  }

  # Add obsp layers
  if (is.character(obsp_layer_names)) {
    obsp <- SOMACollectionCreate(file.path(ms_rna$uri, "obsp"))
    for (layer in obsp_layer_names) {
      obsp$add_new_sparse_ndarray(
        key = layer,
        type = arrow::int32(),
        shape = c(n_obs, n_obs)
      )
      obsp$get(layer)$write(create_sparse_matrix_with_int_dims(
        nrows = n_obs,
        ncols = n_obs
      ))
    }
    obsp$close()
    ms_rna$add_new_collection(obsp, "obsp")
  }

  # Add varp layers
  if (is.character(varp_layer_names)) {
    varp <- SOMACollectionCreate(file.path(ms_rna$uri, "varp"))
    for (layer in varp_layer_names) {
      varp$add_new_sparse_ndarray(
        key = layer,
        type = arrow::int32(),
        shape = c(n_var, n_var)
      )
      varp$get(layer)$write(create_sparse_matrix_with_int_dims(
        nrows = n_var,
        ncols = n_var
      ))
    }
    varp$close()
    ms_rna$add_new_collection(varp, "varp")
  }

  ms_rna$close()

  experiment$ms$set(ms_rna, name = "RNA")
  experiment$ms$close()

  if (is.null(mode)) {
    experiment$close()
  } else if (mode == "READ") {
    experiment$close()
    experiment <- SOMAExperimentOpen(uri, mode = mode)
  }
  experiment
}

# Creates a SOMASparseNDArray with domains of `[0, 2^31 - 1]` and non-zero
# values at `(0,0)`, `(2^31 - 2, 2^31 - 2)` and `(2^31 - 1, 2^31 - 1)`. This is
# intended to test R's ability to read from arrays created with tiledbsoma-py
# before the default domain was changed to `[0, 2^31)`.
#
#  Row/Column:   0      ...    2147483646 |    2147483647
#  0             1      ...    0          |    0
#  ...           ...    ...    ...        |    ...
#  2147483646    0      ...    2          |    0
#  ---------------------------------------|---------------
#  2147483647    0      ...    0          |    3

create_and_populate_32bit_sparse_nd_array <- function(uri) {

  df <- data.frame(
    soma_dim_0 = bit64::as.integer64(c(0, 2^31 - 2, 2^31 - 1)),
    soma_dim_1 = bit64::as.integer64(c(0, 2^31 - 2, 2^31 - 1)),
    soma_data = c(1L, 2L, 3L)
  )

  tdb_dims <- mapply(
    tiledb::tiledb_dim,
    name = c("soma_dim_0", "soma_dim_1"),
    MoreArgs = list(
      domain = c(bit64::as.integer64(0), 2^31 - 1),
      tile = bit64::as.integer64(2),
      type = "INT64"
    ),
    SIMPLIFY = FALSE
  )

  tdb_attr <- tiledb::tiledb_attr(
    name = "soma_data",
    type = "INT32",
  )

  tdb_schema <- tiledb::tiledb_array_schema(
    domain = tiledb::tiledb_domain(tdb_dims),
    attrs = tdb_attr,
    sparse = TRUE
  )

  tiledb::tiledb_array_create(uri, schema = tdb_schema)
  arr <- tiledb::tiledb_array(uri, "WRITE")
  arr[] <- df
  uri
}
