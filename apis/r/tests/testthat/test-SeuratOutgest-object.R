test_that("Load Seurat object from ExperimentQuery mechanics", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("seurat-experiment-query-whole")
  n_obs <- 20L
  n_var <- 10L
  n_pcs <- 50L
  n_umaps <- 2L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    obsm_layers = c(X_pca = n_pcs, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs),
    obsp_layer_names = 'connectivities',
    # No varp in Seurat
    mode = "READ"
  )
  on.exit(experiment$close())

  # Create the query
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )
  # expect_no_condition(obj <- query$to_seurat())
  obj <- query$to_seurat()
  expect_s4_class(obj, 'Seurat')
  expect_identical(dim(obj), c(n_var, n_obs))
  expect_identical(rownames(obj), paste0('feature', query$var_joinids()$as_vector()))
  expect_identical(colnames(obj), paste0('cell', query$obs_joinids()$as_vector()))
  expect_true(all(query$obs_df$attrnames() %in% names(obj[[]])))
  expect_identical(SeuratObject::Assays(obj), 'RNA')
  expect_s4_class(rna <- obj[['RNA']], 'Assay')
  expect_identical(rownames(rna), rownames(obj))
  expect_identical(colnames(rna), colnames(obj))
  expect_identical(
    lapply(list(SeuratObject::Reductions(obj)), sort),
    lapply(list(c('pca', 'umap')), sort)
  )
  expect_s4_class(pca <- obj[['pca']], 'DimReduc')
  expect_identical(SeuratObject::Cells(pca), colnames(obj))
  expect_identical(rownames(SeuratObject::Loadings(pca)), rownames(obj))
  expect_identical(ncol(pca), n_pcs)
  expect_s4_class(umap <- obj[['umap']], 'DimReduc')
  expect_identical(SeuratObject::Cells(umap), colnames(obj))
  expect_identical(ncol(umap), n_umaps)
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(umap)))
  expect_identical(SeuratObject::Graphs(obj), 'connectivities')
  expect_s4_class(graph <- obj[['connectivities']], 'Graph')
  expect_identical(dim(graph), c(n_obs, n_obs))
  expect_identical(rownames(graph), colnames(obj))
  expect_identical(colnames(graph), colnames(obj))

  # Test named
  # expect_no_condition(obj <- query$to_seurat(
  #   obs_index = 'baz',
  #   var_index = 'quux'
  # ))
  obj <- query$to_seurat(
    obs_index = 'baz',
    var_index = 'quux'
  )
  expect_s4_class(obj, 'Seurat')
  expect_identical(dim(obj), c(n_var, n_obs))
  expect_identical(
    rownames(obj),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_identical(
    colnames(obj),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::Assays(obj), 'RNA')
  expect_false(all(query$obs_df$attrnames() %in% names(obj[[]])))
  expect_true(all(setdiff(query$obs_df$attrnames(), 'baz') %in% names(obj[[]])))
  expect_s4_class(rna <- obj[['RNA']], 'Assay')
  expect_identical(rownames(rna), rownames(obj))
  expect_identical(colnames(rna), colnames(obj))
  expect_identical(
    lapply(list(SeuratObject::Reductions(obj)), sort),
    lapply(list(c('pca', 'umap')), sort)
  )
  expect_s4_class(pca <- obj[['pca']], 'DimReduc')
  expect_identical(SeuratObject::Cells(pca), colnames(obj))
  expect_identical(rownames(SeuratObject::Loadings(pca)), rownames(obj))
  expect_identical(ncol(pca), n_pcs)
  expect_s4_class(umap <- obj[['umap']], 'DimReduc')
  expect_identical(SeuratObject::Cells(umap), colnames(obj))
  expect_identical(ncol(umap), n_umaps)
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(umap)))
  expect_identical(SeuratObject::Graphs(obj), 'connectivities')
  expect_s4_class(graph <- obj[['connectivities']], 'Graph')
  expect_identical(dim(graph), c(n_obs, n_obs))
  expect_identical(rownames(graph), colnames(obj))
  expect_identical(colnames(graph), colnames(obj))

  # Test `X_layers`
  # expect_no_condition(obj <- query$to_seurat(c(counts = 'counts')))
  obj <- query$to_seurat(c(counts = 'counts'))
  expect_s4_class(
    counts <- SeuratObject::GetAssayData(obj[['RNA']], 'counts'),
    'dgCMatrix'
  )
  expect_s4_class(
    data <- SeuratObject::GetAssayData(obj[['RNA']], 'data'),
    'dgCMatrix'
  )
  expect_identical(counts, data)
  # expect_no_condition(obj <- query$to_seurat(c(data = 'logcounts')))
  obj <- query$to_seurat(c(data = 'logcounts'))
  expect_s4_class(
    SeuratObject::GetAssayData(obj[['RNA']], 'data'),
    'dgCMatrix'
  )
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::GetAssayData(
    obj[['RNA']],
    'counts'
  )))

  # Test suppress reductions
  # expect_no_condition(obj <- query$to_seurat(obsm_layers = FALSE))
  obj <- query$to_seurat(obsm_layers = FALSE)
  expect_length(SeuratObject::Reductions(obj), 0L)
  # expect_no_condition(obj <- query$to_seurat(obsm_layers = NA))
  obj <- query$to_seurat(obsm_layers = NA)
  expect_length(SeuratObject::Reductions(obj), 0L)
  # expect_no_condition(obj <- query$to_seurat(obsm_layers = 'umap'))
  obj <- query$to_seurat(obsm_layers = 'umap')
  expect_identical(SeuratObject::Reductions(obj), 'umap')
  expect_error(obj[['pca']])

  # Test suppress loadings
  # expect_no_condition(obj <- query$to_seurat(varm_layers = FALSE))
  obj <- query$to_seurat(varm_layers = FALSE)
  expect_identical(
    lapply(list(SeuratObject::Reductions(obj)), sort),
    lapply(list(c('pca', 'umap')), sort)
  )
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(obj[['pca']])))

  # Test suppress graphs
  # expect_no_condition(obj <- query$to_seurat(obsp_layers = FALSE))
  obj <- query$to_seurat(obsp_layers = FALSE)
  expect_length(SeuratObject::Graphs(obj), 0L)

  # Test suppress cell-level meta data
  # expect_no_condition(obj <- query$to_seurat(obs_column_names = FALSE))
  obj <- query$to_seurat(obs_column_names = FALSE)
  expect_false(any(query$obs_df$attrnames() %in% names(obj[[]])))

  # Test `X_layers` assertions
  expect_error(query$to_seurat(NULL))
  expect_error(query$to_seurat(FALSE))
  expect_error(query$to_seurat(1))
  expect_error(query$to_seurat('counts'))
  expect_error(query$to_seurat(unlist(list(counts = 'counts', 'logcounts'))))
  expect_error(query$to_seurat(list(counts = 'counts', data = 'logcounts')))
  expect_error(query$to_seurat(c(a = 'counts')))
  expect_error(query$to_seurat(c(scale.data = 'counts')))
  expect_error(query$to_seurat(c(data = 'tomato')))

  # Test `obs_index` assertions
  expect_error(query$to_seurat(obs_index = FALSE))
  expect_error(query$to_seurat(obs_index = NA_character_))
  expect_error(query$to_seurat(obs_index = 1))
  expect_error(query$to_seurat(obs_index = c('baz', 'foo')))
  expect_error(query$to_seurat(obs_index = 'tomato'))

  # Test `obs_column_names` assertions
  expect_error(query$to_seurat(obs_column_names = 1L))
  expect_error(query$to_seurat(obs_column_names = c(
    NA_character_,
    NA_character_
  )))
  expect_error(query$to_seurat(obs_column_names = c(TRUE, FALSE)))
  expect_error(query$to_seurat(obs_column_names = 'tomato'))

  # Test `obsm_layers` assertions
  expect_error(query$to_seurat(obsm_layers = 1L))
  expect_warning(query$to_seurat(obsm_layers = 'tomato'))

  # Test `varm_layers` assertions
  expect_error(query$to_seurat(varm_layers = 1L))
  expect_error(query$to_seurat(varm_layers = 'PCs'))
  expect_warning(query$to_seurat(varm_layers = c(tomato = 'PCs')))
  expect_warning(query$to_seurat(varm_layers = c(X_pca = 'tomato')))

  # Test `obsp_layers` assertions
  expect_error(query$to_seurat(obsp_layers = 1L))
  expect_warning(query$to_seurat(obsp_layers = 'tomato'))
})

test_that("Load Seurat object from sliced ExperimentQuery", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("seurat-experiment-query-sliced")
  n_obs <- 1001L
  n_var <- 99L
  n_pcs <- 50L
  n_umaps <- 2L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    obsm_layers = c(X_pca = n_pcs, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs),
    obsp_layer_names = 'connectivities',
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
  n_var_slice <- length(var_slice)
  n_obs_slice <- length(obs_slice)
  # expect_no_condition(obj <- query$to_seurat())
  obj <- query$to_seurat()
  expect_s4_class(obj, 'Seurat')
  expect_identical(dim(obj), c(n_var_slice, n_obs_slice))
  expect_identical(rownames(obj), paste0('feature', query$var_joinids()$as_vector()))
  expect_identical(colnames(obj), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(
    lapply(list(names(obj)), sort),
    lapply(list(c('RNA', 'connectivities', 'pca', 'umap')), sort)
  )

  # Test named
  # expect_no_condition(obj <- query$to_seurat(
  #   obs_index = 'baz',
  #   var_index = 'quux'
  # ))
  obj <- query$to_seurat(
    obs_index = 'baz',
    var_index = 'quux'
  )
  expect_s4_class(obj, 'Seurat')
  expect_identical(dim(obj), c(n_var_slice, n_obs_slice))
  expect_identical(
    rownames(obj),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_identical(
    colnames(obj),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    lapply(list(names(obj)), sort),
    lapply(list(c('RNA', 'connectivities', 'pca', 'umap')), sort)
  )
})

test_that("Load Seurat object from indexed ExperimentQuery", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("seurat-experiment-query-value-filters")
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
    X_layer_names = c("counts", "logcounts"),
    obsm_layers = c(X_pca = n_pcs, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs),
    obsp_layer_names = 'connectivities',
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
  n_var_select <- length(var_label_values)
  n_obs_select <- length(obs_label_values)
  # expect_no_condition(obj <- query$to_seurat())
  obj <- query$to_seurat()
  expect_s4_class(obj, 'Seurat')
  expect_identical(dim(obj), c(n_var_select, n_obs_select))
  expect_identical(rownames(obj), paste0('feature', query$var_joinids()$as_vector()))
  expect_identical(colnames(obj), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(
    lapply(list(names(obj)), sort),
    lapply(list(c('RNA', 'connectivities', 'pca', 'umap')), sort)
  )

  # Test named
  # expect_no_condition(obj <- query$to_seurat(
  #   obs_index = 'baz',
  #   var_index = 'quux'
  # ))
  obj <- query$to_seurat(
    obs_index = 'baz',
    var_index = 'quux'
  )
  expect_s4_class(obj, 'Seurat')
  expect_identical(dim(obj), c(n_var_select, n_obs_select))
  expect_identical(
    rownames(obj),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_identical(
    colnames(obj),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    lapply(list(names(obj)), sort),
    lapply(list(c('RNA', 'connectivities', 'pca', 'umap')), sort)
  )
})
