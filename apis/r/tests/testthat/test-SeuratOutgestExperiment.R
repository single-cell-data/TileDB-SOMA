test_that("Load Seurat object from SOMA Experiment mechanics", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("seurat-experiment-whole")
  n_obs <- 1001L
  n_var <- 99L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )
  experiment <- SOMAExperimentOpen(experiment$uri, mode = 'WRITE')
  exp_ms_rna <- experiment$ms$get('RNA')
  # Add embeddings
  n_pcs <- 50L
  n_umaps <- 2L
  obsm <- SOMACollectionCreate(file.path(exp_ms_rna$uri, 'obsm'))
  obsm$add_new_sparse_ndarray(
    key = 'X_pca',
    type = arrow::int32(),
    shape = c(n_obs, n_pcs)
  )
  obsm$get('X_pca')$write(create_sparse_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_pcs
  ))
  obsm$add_new_sparse_ndarray(
    key = 'X_umap',
    type = arrow::int32(),
    shape = c(n_obs, n_umaps)
  )
  obsm$get('X_umap')$write(create_sparse_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_umaps,
    seed = 2L
  ))
  obsm$close()
  exp_ms_rna$add_new_collection(obsm, 'obsm')
  # Add loadings
  varm <- SOMACollectionCreate(file.path(exp_ms_rna$uri, 'varm'))
  varm$add_new_sparse_ndarray(
    key = 'PCs',
    type = arrow::int32(),
    shape = c(n_var, n_pcs)
  )
  varm$get('PCs')$write(create_sparse_matrix_with_int_dims(
    nrows = n_var,
    ncols = n_pcs
  ))
  varm$close()
  exp_ms_rna$add_new_collection(varm, 'varm')
  # Add graph
  obsp <- SOMACollectionCreate(file.path(exp_ms_rna$uri, 'obsp'))
  obsp$add_new_sparse_ndarray(
    key = 'connectivities',
    type = arrow::int32(),
    shape = c(n_obs, n_obs)
  )
  obsp$get('connectivities')$write(create_sparse_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_obs
  ))
  obsp$close()
  exp_ms_rna$add_new_collection(obsp, 'obsp')
  # Close and reopen in read mode
  exp_ms_rna$close()
  experiment$close()
  experiment <- SOMAExperimentOpen(experiment$uri)
  var_tbl <- experiment$ms$get('RNA')$get('var')$read()
  obs_tbl <- experiment$obs$read()
  X_layers <- list(RNA = c(counts = 'counts', data = 'logcounts'))
  expect_no_condition(obj <- experiment$to_seurat(X_layers))
  expect_s4_class(obj, 'Seurat')
  expect_identical(dim(obj), c(n_var, n_obs))
  expect_identical(
    rownames(obj),
    paste0(
      'feature',
      var_tbl$GetColumnByName('soma_joinid')$as_vector()
    )
  )
  expect_identical(
    colnames(obj),
    paste0(
      'cell',
      obs_tbl$GetColumnByName('soma_joinid')$as_vector()
    )
  )
  expect_true(all(experiment$obs$attrnames() %in% names(obj[[]])))
  expect_identical(SeuratObject::Assays(obj), 'RNA')
  expect_s4_class(rna <- obj[['RNA']], 'Assay')
  expect_identical(rownames(rna), rownames(obj))
  expect_identical(colnames(rna), colnames(obj))
  expect_identical(sort(SeuratObject::Reductions(obj)), c('pca', 'umap'))
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
  expect_no_condition(obj <- experiment$to_seurat(
    X_layers,
    obs_index = 'baz',
    var_index = c(RNA = 'quux')
  ))
  expect_s4_class(obj, 'Seurat')
  expect_identical(dim(obj), c(n_var, n_obs))
  expect_identical(
    rownames(obj),
    var_tbl$GetColumnByName('quux')$as_vector()
  )
  expect_identical(
    colnames(obj),
    obs_tbl$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::Assays(obj), 'RNA')
  expect_false(all(experiment$obs$attrnames() %in% names(obj[[]])))
  expect_true(all(setdiff(experiment$obs$attrnames(), 'baz') %in% names(obj[[]])))
  expect_s4_class(rna <- obj[['RNA']], 'Assay')
  expect_identical(rownames(rna), rownames(obj))
  expect_identical(colnames(rna), colnames(obj))
  expect_identical(sort(SeuratObject::Reductions(obj)), c('pca', 'umap'))
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
  expect_no_condition(obj <- experiment$to_seurat(list(RNA = c(counts = 'counts'))))
  expect_s4_class(
    counts <- SeuratObject::GetAssayData(obj[['RNA']], 'counts'),
    'dgCMatrix'
  )
  expect_s4_class(
    data <- SeuratObject::GetAssayData(obj[['RNA']], 'data'),
    'dgCMatrix'
  )
  expect_identical(counts, data)
  expect_no_condition(obj <- experiment$to_seurat(
    list(RNA = c(data = 'logcounts'))
  ))
  expect_s4_class(
    SeuratObject::GetAssayData(obj[['RNA']], 'data'),
    'dgCMatrix'
  )
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::GetAssayData(
    obj[['RNA']],
    'counts'
  )))
  # Test suppress reductions
  expect_no_condition(obj <- experiment$to_seurat(X_layers, obsm_layers = FALSE))
  expect_length(SeuratObject::Reductions(obj), 0L)
  expect_no_condition(obj <- experiment$to_seurat(X_layers, obsm_layers = NA))
  expect_length(SeuratObject::Reductions(obj), 0L)
  expect_no_condition(obj <- experiment$to_seurat(
    X_layers,
    obsm_layers = list(RNA = 'umap')
  ))
  expect_identical(SeuratObject::Reductions(obj), 'umap')
  expect_error(obj[['pca']])
  # Test suppress loadings
  expect_no_condition(obj <- experiment$to_seurat(X_layers, varm_layers = FALSE))
  expect_identical(sort(SeuratObject::Reductions(obj)), c('pca', 'umap'))
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(obj[['pca']])))
  # Test suppress graphs
  expect_no_condition(obj <- experiment$to_seurat(X_layers, obsp_layers = FALSE))
  expect_length(SeuratObject::Graphs(obj), 0L)
  # Test suppress cell-level meta data
  expect_no_condition(obj <- experiment$to_seurat(
    X_layers,
    obs_column_names = FALSE
  ))
  expect_false(any(experiment$obs$attrnames() %in% names(obj[[]])))
  # Test `X_layers` assertions
  expect_error(experiment$to_seurat(NULL))
  expect_error(experiment$to_seurat(FALSE))
  expect_error(experiment$to_seurat(1))
  expect_error(experiment$to_seurat('counts'))
  expect_error(experiment$to_seurat(unlist(list(counts = 'counts', 'logcounts'))))
  expect_error(experiment$to_seurat(list(counts = 'counts', data = 'logcounts')))
  expect_error(experiment$to_seurat(c(a = 'counts')))
  expect_error(experiment$to_seurat(c(scale.data = 'counts')))
  expect_error(query$to_seurat(c(data = 'tomato')))
  # Test `obs_index` assertions
  expect_error(experiment$to_seurat(X_layers, obs_index = FALSE))
  expect_error(experiment$to_seurat(X_layers, obs_index = NA_character_))
  expect_error(experiment$to_seurat(X_layers, obs_index = 1))
  expect_error(experiment$to_seurat(X_layers, obs_index = c('baz', 'foo')))
  expect_error(experiment$to_seurat(X_layers, obs_index = 'tomato'))
  # Test `obs_column_names` assertions
  expect_error(experiment$to_seurat(X_layers, obs_column_names = 1L))
  expect_error(experiment$to_seurat(
    X_layers,
    obs_column_names = c(
      NA_character_,
      NA_character_
    )
  ))
  expect_error(experiment$to_seurat(X_layers, obs_column_names = c(TRUE, FALSE)))
  expect_error(experiment$to_seurat(X_layers, obs_column_names = 'tomato'))
  expect_error(experiment$to_seurat(X_layers, obs_column_names = 'tomato'))
  # Test `obsm_layers` assertions
  expect_error(experiment$to_seurat(X_layers, obsm_layers = 1L))
  expect_error(experiment$to_seurat(X_layers, obsm_layers = 'tomato'))
  expect_error(experiment$to_seurat(X_layers, obsm_layers = 'umap'))
  # Test `varm_layers` assertions
  expect_error(experiment$to_seurat(X_layers, varm_layers = 1L))
  expect_error(experiment$to_seurat(X_layers, varm_layers = 'PCs'))
  expect_error(experiment$to_seurat(X_layers, varm_layers = c(tomato = 'PCs')))
  expect_error(experiment$to_seurat(X_layers, varm_layers = c(X_pca = 'tomato')))
  # Test `obsp_layers` assertions
  expect_error(experiment$to_seurat(X_layers, obsp_layers = 1L))
  expect_error(experiment$to_seurat(X_layers, obsp_layers = 'tomato'))
})
