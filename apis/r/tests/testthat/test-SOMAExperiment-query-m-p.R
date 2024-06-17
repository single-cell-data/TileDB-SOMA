test_that("Load *m and *p layers from SOMAExperimentAxisQuery mechanics", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern="m-p-experiment-query-mechanics")

  n_obs <- 20L
  n_var <- 10L
  n_pcs <- 50L
  n_ics <- 30L
  n_umaps <- 2L

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = "counts",
    obsm_layers = c(X_pca = n_pcs, 'dense:X_ica' = n_ics, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs, "dense:ICs" = n_ics),
    obsp_layer_names = 'connectivities',
    varp_layer_names = 'network',
    mode = 'READ'
  )
  on.exit(experiment$close(), add = TRUE)

  # Create the query
  query <- SOMAExperimentAxisQuery$new(experiment, measurement_name = "RNA")

  # Read `obsm` layers
  expect_s3_class(pc_qry <- query$obsm("X_pca"), "SOMASparseNDArrayRead")
  pc_mat <- pc_qry$sparse_matrix()$concat()
  expect_s4_class(pc_mat, "dgTMatrix")
  expect_identical(dim(pc_mat), c(n_obs, n_pcs))

  expect_s3_class(umap_qry <- query$obsm("X_umap"), "SOMASparseNDArrayRead")
  umap_mat <- umap_qry$sparse_matrix()$concat()
  expect_s4_class(umap_mat, "dgTMatrix")
  expect_identical(dim(umap_mat), c(n_obs, n_umaps))

  expect_error(query$obsm("X_ica"), regexp = "obsm layer must be a SOMASparseNDArray")

  # Read `varm` layers
  expect_s3_class(pc_qry <- query$varm("PCs"), "SOMASparseNDArrayRead")
  pc_mat <- pc_qry$sparse_matrix()$concat()
  expect_s4_class(pc_mat, "dgTMatrix")
  expect_identical(dim(pc_mat), c(n_var, n_pcs))

  expect_error(query$varm("ICs"), regexp = "varm layer must be a SOMASparseNDArray")

  # Read `obsp` layers
  expect_s3_class(con_qry <- query$obsp("connectivities"), "SOMASparseNDArrayRead")
  con_mat <- con_qry$sparse_matrix()$concat()
  expect_s4_class(con_mat, "dgTMatrix")
  expect_identical(dim(con_mat), c(n_obs, n_obs))

  # Read `varp` layers
  expect_s3_class(net_qry <- query$varp("network"), "SOMASparseNDArrayRead")
  net_mat <- net_qry$sparse_matrix()$concat()
  expect_s4_class(net_mat, "dgTMatrix")
  expect_identical(dim(net_mat), c(n_var, n_var))

  # Test assertions
  expect_error(query$obsm())
  expect_error(query$obsm(1L))
  expect_error(query$obsm(c("X_pca", "X_umap")))
  expect_error(query$obsm("missing-layer"))

  expect_error(query$varm())
  expect_error(query$varm(1L))
  expect_error(query$varm(c("X_pca", "X_umap")))
  expect_error(query$varm("missing-layer"))

  expect_error(query$obsp())
  expect_error(query$obsp(1L))
  expect_error(query$obsp(c("X_pca", "X_umap")))
  expect_error(query$obsp("missing-layer"))

  expect_error(query$varp())
  expect_error(query$varp(1L))
  expect_error(query$varp(c("X_pca", "X_umap")))
  expect_error(query$varp("missing-layer"))
})

test_that("SOMAExperimentAxisQuery without *m and *p layers mechanics", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern="m-p-missing-experiment-query-mechanics")

  n_obs <- 1001L
  n_var <- 99L

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = "counts",
    mode = 'READ'
  )
  on.exit(experiment$close(), add = TRUE)

  # Create the query
  query <- SOMAExperimentAxisQuery$new(experiment, measurement_name = "RNA")

  expect_error(query$obsm(), class = "noLayersError")
  expect_error(query$varm(), class = "noLayersError")
  expect_error(query$obsp(), class = "noLayersError")
  expect_error(query$varp(), class = "noLayersError")
})

test_that("Load *m and *p layers from sliced SOMAExperimentAxisQuery", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern="m-p-experiment-query-sliced")

  n_obs <- 1001L
  n_var <- 99L
  n_pcs <- 50L
  n_umaps <- 2L

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = "counts",
    obsm_layers = c(X_pca = n_pcs, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs),
    obsp_layer_names = 'connectivities',
    varp_layer_names = 'network',
    mode = 'READ'
  )
  on.exit(experiment$close(), add = TRUE)

  # Create the query
  obs_slice <- bit64::as.integer64(seq(3, 72))
  var_slice <- bit64::as.integer64(seq(7, 21))
  n_var_slice <- length(var_slice)
  n_obs_slice <- length(obs_slice)
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(coords = list(soma_joinid = obs_slice)),
    var_query = SOMAAxisQuery$new(coords = list(soma_joinid = var_slice))
  )

  # Read `obsm` layers
  expect_s3_class(pc_qry <- query$obsm("X_pca"), "SOMASparseNDArrayRead")
  pc_tbl <- pc_qry$tables()$concat()
  expect_s3_class(pc_tbl, "Table")
  pc_obs <- pc_tbl$GetColumnByName("soma_dim_0")$as_vector()
  expect_gte(min(pc_obs), min(obs_slice))
  expect_lte(max(pc_obs), max(obs_slice))

  expect_s3_class(umap_qry <- query$obsm("X_umap"), "SOMASparseNDArrayRead")
  umap_tbl <- umap_qry$tables()$concat()
  expect_s3_class(umap_tbl, "Table")
  umap_obs <- umap_tbl$GetColumnByName("soma_dim_0")$as_vector()
  expect_gte(min(umap_obs), min(obs_slice))
  expect_lte(max(umap_obs), max(obs_slice))

  # Read `varm` layers
  expect_s3_class(pc_qry <- query$varm("PCs"), "SOMASparseNDArrayRead")
  pc_tbl <- pc_qry$tables()$concat()
  expect_s3_class(pc_tbl, "Table")
  pc_var <- pc_tbl$GetColumnByName("soma_dim_0")$as_vector()
  expect_gte(min(pc_var), min(var_slice))
  expect_lte(max(pc_var), max(var_slice))

  # Read `obsp` layers
  expect_s3_class(con_qry <- query$obsp("connectivities"), "SOMASparseNDArrayRead")
  con_tbl <- con_qry$tables()$concat()
  expect_s3_class(con_tbl, "Table")
  for (axis in paste0("soma_dim_", 0:1)) {
    con_dim <- con_tbl$GetColumnByName(axis)$as_vector()
    expect_gte(
      min(con_dim),
      min(obs_slice),
      label = axis,
      expected.label = 'lower obsp range'
    )
    expect_lte(
      max(con_dim),
      max(obs_slice),
      label = axis,
      expected.label = 'upper obsp range'
    )
  }

  # Read `varp` layers
  expect_s3_class(net_qry <- query$varp("network"), "SOMASparseNDArrayRead")
  net_tbl <- net_qry$tables()$concat()
  expect_s3_class(net_tbl, "Table")
  for (axis in paste0("soma_dim_", 0:1)) {
    net_dim <- net_tbl$GetColumnByName(axis)$as_vector()
    expect_gte(
      min(net_dim),
      min(var_slice),
      label = axis,
      expected.label = 'lower varp range'
    )
    expect_lte(
      max(net_dim),
      max(var_slice),
      label = axis,
      expected.label = 'upper varp range'
    )
  }

})

test_that("Load *m and *p layers from indexed SOMAExperimentAxisQuery", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern="m-p-experiment-query-indexed")

  n_obs <- 1001L
  n_var <- 99L
  n_pcs <- 50L
  n_umaps <- 2L
  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")

  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = "counts",
    obsm_layers = c(X_pca = n_pcs, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs),
    obsp_layer_names = 'connectivities',
    varp_layer_names = 'network',
    mode = 'READ'
  )
  on.exit(experiment$close(), add = TRUE)

  # Create the query
  obs_value_filter <- paste0(
    sprintf("baz == '%s'", obs_label_values),
    collapse = "||"
  )
  var_value_filter <- paste0(
    sprintf("quux == '%s'", var_label_values),
    collapse = "||"
  )
  n_obs_select <- length(obs_label_values)
  n_var_select <- length(var_label_values)
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(value_filter = obs_value_filter),
    var_query = SOMAAxisQuery$new(value_filter = var_value_filter)
  )


  obs_ids <- query$obs("soma_joinid")$concat()$GetColumnByName("soma_joinid")$as_vector()
  var_ids <- query$var("soma_joinid")$concat()$GetColumnByName("soma_joinid")$as_vector()
  expect_length(obs_ids, n_obs_select)
  expect_length(var_ids, n_var_select)

  # Read `obsm` layers
  expect_s3_class(pc_qry <- query$obsm("X_pca"), "SOMASparseNDArrayRead")
  pc_tbl <- pc_qry$tables()$concat()
  expect_s3_class(pc_tbl, "Table")
  pc_obs <- pc_tbl$GetColumnByName("soma_dim_0")$as_vector()
  expect_in(pc_obs, obs_ids)

  expect_s3_class(umap_qry <- query$obsm("X_umap"), "SOMASparseNDArrayRead")
  umap_tbl <- umap_qry$tables()$concat()
  expect_s3_class(umap_tbl, "Table")
  umap_obs <- umap_tbl$GetColumnByName("soma_dim_0")$as_vector()
  expect_in(umap_obs, obs_ids)

  # Read `varm` layers
  expect_s3_class(pc_qry <- query$varm("PCs"), "SOMASparseNDArrayRead")
  pc_tbl <- pc_qry$tables()$concat()
  expect_s3_class(pc_tbl, "Table")
  pc_var <- pc_tbl$GetColumnByName("soma_dim_0")$as_vector()
  expect_in(pc_var, var_ids)

  # Read `obsp` layers
  expect_s3_class(con_qry <- query$obsp("connectivities"), "SOMASparseNDArrayRead")
  con_tbl <- con_qry$tables()$concat()
  expect_s3_class(con_tbl, "Table")
  for (axis in paste0("soma_dim_", 0:1)) {
    con_dim <- con_tbl$GetColumnByName(axis)$as_vector()
    expect_in(con_dim, obs_ids)
  }

  # Read `varp` layers
  expect_s3_class(net_qry <- query$varp("network"), "SOMASparseNDArrayRead")
  net_tbl <- net_qry$tables()$concat()
  expect_s3_class(net_tbl, "Table")
  for (axis in paste0("soma_dim_", 0:1)) {
    net_dim <- net_tbl$GetColumnByName(axis)$as_vector()
    expect_in(net_dim, var_ids)
  }

})
