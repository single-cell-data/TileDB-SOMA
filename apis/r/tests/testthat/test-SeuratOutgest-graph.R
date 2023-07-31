test_that("Load graph from ExperimentQuery mechanics", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("graph-experiment-query-whole")
  n_obs <- 20L
  n_var <- 10L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    obsp_layer_names = "connectivities",
    # No varp in Seurat
    mode = "READ"
  )
  on.exit(experiment$close())

  # Create the query
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )
  expect_no_condition(graph <- query$to_seurat_graph('connectivities'))
  expect_s4_class(graph, 'Graph')
  expect_identical(dim(graph), c(n_obs, n_obs))
  expect_identical(rownames(graph), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(colnames(graph), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(SeuratObject::DefaultAssay(graph), 'RNA')

  # Test named
  expect_no_condition(named <- query$to_seurat_graph('connectivities', 'baz'))
  expect_s4_class(named, 'Graph')
  expect_identical(dim(named), c(n_obs, n_obs))
  expect_identical(
    rownames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    colnames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::DefaultAssay(named), 'RNA')

  # Test `graph` assertions
  expect_error(query$to_seurat_graph(NULL))
  expect_error(query$to_seurat_graph(TRUE))
  expect_error(query$to_seurat_graph(NA))
  expect_error(query$to_seurat_graph(1))
  expect_error(query$to_seurat_graph(c('connectivities',)))
  expect_error(query$to_seurat_graph('tomato'))

  # Test `obs_index` assertions
  expect_error(query$to_seurat_graph('connectivities', obs_index = FALSE))
  expect_error(query$to_seurat_graph(
    'connectivities',
    obs_index = NA_character_)
  )
  expect_error(query$to_seurat_graph('connectivities', obs_index = 1))
  expect_error(query$to_seurat_graph(
    'connectivities',
    obs_index = c('baz', 'foo'))
  )
  expect_error(query$to_seurat_graph('connectivities', obs_index = 'tomato'))
})

test_that("Load graph from sliced ExperimentQuery", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("graph-experiment-query-sliced")
  n_obs <- 1001L
  n_var <- 99L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    obsp_layer_names = "connectivities",
    # No varp in Seurat
    mode = "READ"
  )
  on.exit(experiment$close())

  # Create the query
  obs_slice <- bit64::as.integer64(seq(3, 72))
  var_slice <- bit64::as.integer64(seq(7, 21))
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(coords = list(soma_joinid = obs_slice)),
    var_query = SOMAAxisQuery$new(coords = list(soma_joinid = var_slice))
  )
  n_obs_slice <- length(obs_slice)
  expect_no_condition(graph <- query$to_seurat_graph('connectivities'))
  expect_s4_class(graph, 'Graph')
  expect_identical(dim(graph), c(n_obs_slice, n_obs_slice))
  expect_identical(rownames(graph), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(colnames(graph), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(SeuratObject::DefaultAssay(graph), 'RNA')

  # Test named
  expect_no_condition(named <- query$to_seurat_graph('connectivities', 'baz'))
  expect_s4_class(named, 'Graph')
  expect_identical(dim(named), c(n_obs_slice, n_obs_slice))
  expect_identical(
    rownames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    colnames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::DefaultAssay(named), 'RNA')
})

test_that("Load graph from indexed ExperimentQuery", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("graph-experiment-query-value-filters")
  n_obs <- 1001L
  n_var <- 99L
  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    obsp_layer_names = "connectivities",
    # No varp in Seurat
    mode = "READ"
  )
  on.exit(experiment$close())

  # Create the query
  obs_value_filter <- paste0(
    sprintf("baz == '%s'", obs_label_values),
    collapse = "||"
  )
  var_value_filter <- paste0(
    sprintf("quux == '%s'", var_label_values),
    collapse = "||"
  )
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(value_filter = obs_value_filter),
    var_query = SOMAAxisQuery$new(value_filter = var_value_filter)
  )
  n_obs_select <- length(obs_label_values)
  expect_no_condition(graph <- query$to_seurat_graph('connectivities'))
  expect_s4_class(graph, 'Graph')
  expect_identical(dim(graph), c(n_obs_select, n_obs_select))
  expect_identical(rownames(graph), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(colnames(graph), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(SeuratObject::DefaultAssay(graph), 'RNA')

  # Test named
  expect_no_condition(named <- query$to_seurat_graph('connectivities', 'baz'))
  expect_s4_class(named, 'Graph')
  expect_identical(dim(named), c(n_obs_select, n_obs_select))
  expect_identical(
    rownames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(rownames(named), obs_label_values)
  expect_identical(
    colnames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(colnames(named), obs_label_values)
  expect_identical(SeuratObject::DefaultAssay(named), 'RNA')
})
