test_that("Load assay from ExperimentQuery mechanics", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("assay-experiment-query-whole")
  n_obs <- 20L
  n_var <- 10L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )
  expect_no_condition(assay <- query$to_seurat_assay())
  expect_s4_class(assay, 'Assay')
  expect_identical(dim(assay), c(n_var, n_obs))
  expect_s4_class(SeuratObject::GetAssayData(assay, 'counts'), 'dgCMatrix')
  expect_s4_class(SeuratObject::GetAssayData(assay, 'data'), 'dgCMatrix')
  scale.data <- SeuratObject::GetAssayData(assay, 'scale.data')
  expect_true(is.matrix(scale.data))
  expect_equal(dim(scale.data), c(0, 0))
  expect_equal(SeuratObject::Key(assay), 'rna_')
  expect_equal(names(assay[[]]), query$var_df$attrnames())
  expect_equal(rownames(assay), paste0('feature', seq_len(n_var) - 1L))
  expect_equal(rownames(assay), paste0('feature', query$var_joinids()))
  expect_equal(colnames(assay), paste0('cell', seq_len(n_obs) - 1L))
  expect_equal(colnames(assay), paste0('cell', query$obs_joinids()))
  # Test no counts
  expect_no_condition(nocounts <- query$to_seurat_assay(c(data = 'logcounts')))
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::GetAssayData(
    nocounts,
    'counts'
  )))
  # Test no data (populate `data` with `counts`)
  expect_no_condition(nodata <- query$to_seurat_assay(c(counts = 'counts')))
  expect_identical(
    SeuratObject::GetAssayData(nodata, 'data'),
    SeuratObject::GetAssayData(nodata, 'counts')
  )
  # Test adding `scale.data`
  expect_no_condition(sd <- query$to_seurat_assay(c(
    data = 'logcounts', scale.data = 'counts'
  )))
  expect_s4_class(scaled <- SeuratObject::GetAssayData(sd, 'scale.data'), NA)
  expect_true(is.matrix(scaled))
  expect_equal(dim(scaled), c(n_var, n_obs))
  # Test modifying feature-level meta data
  expect_no_condition(nomf <- query$to_seurat_assay(var_column_names = FALSE))
  expect_equal(dim(nomf[[]]), c(n_var, 0L))
  expect_no_condition(nomf2 <- query$to_seurat_assay(var_column_names = NA))
  expect_equal(dim(nomf2[[]]), c(n_var, 0L))
  # Test using cell and feature names
  expect_no_condition(named <- query$to_seurat_assay(
    cells_index = 'baz',
    features_index = 'quux'
  ))
  expect_identical(
    colnames(named),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(named),
    query$var('quux')$GetColumnByName('quux')$as_vector()
  )
  # Test `X_layers` assertions
  expect_error(query$to_seurat_assay(NULL))
  expect_error(query$to_seurat_assay(FALSE))
  expect_error(query$to_seurat_assay(1))
  expect_error(query$to_seurat_assay('counts'))
  expect_error(query$to_seurat_assay(unlist(list(
    counts = 'counts',
    'logcounts'
  ))))
  expect_error(query$to_seurat_assay(list(
    counts = 'counts',
    data = 'logcounts'
  )))
  expect_error(query$to_seurat_assay(c(a = 'counts')))
  expect_error(query$to_seurat_assay(c(scale.data = 'counts')))
  expect_error(query$to_seurat_assay(c(data = 'tomato')))
  expect_warning(query$to_seurat_assay(c(counts = 'counts', data = 'tomato')))
  # Test `cells_index` assertions
  expect_error(query$to_seurat_assay(cells_index = FALSE))
  expect_error(query$to_seurat_assay(cells_index = NA_character_))
  expect_error(query$to_seurat_assay(cells_index = 1))
  expect_error(query$to_seurat_assay(cells_index = c('baz', 'foo')))
  expect_error(query$to_seurat_assay(cells_index = 'tomato'))
  # Test `features_index` assertions
  expect_error(query$to_seurat_assay(features_index = FALSE))
  expect_error(query$to_seurat_assay(features_index = NA_character_))
  expect_error(query$to_seurat_assay(features_index = 1))
  expect_error(query$to_seurat_assay(features_index = c('baz', 'foo')))
  expect_error(query$to_seurat_assay(features_index = 'tomato'))
  # Test `var_column_names` assertions
  expect_error(query$to_seurat_assay(var_column_names = 1L))
  expect_error(query$to_seurat_assay(var_column_names = c(
    NA_character_,
    NA_character_
  )))
  expect_error(query$to_seurat_assay(var_column_names = c(TRUE, FALSE)))
  expect_error(query$to_seurat_assay(var_column_names = 'tomato'))
})

test_that("Load assay from sliced ExperimentQuery", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("assay-experiment-query-sliced")
  n_obs <- 1001L
  n_var <- 99L
  obs_slice <- bit64::as.integer64(seq(3, 72))
  var_slice <- bit64::as.integer64(seq(7, 21))
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(coords = list(soma_joinid = obs_slice)),
    var_query = SOMAAxisQuery$new(coords = list(soma_joinid = var_slice))
  )
  expect_no_condition(assay <- query$to_seurat_assay())
  expect_s4_class(assay, 'Assay')
  expect_identical(dim(assay), c(length(var_slice), length(obs_slice)))
  expect_identical(rownames(assay), paste0('feature', query$var_joinids()))
  expect_identical(colnames(assay), paste0('cell', query$obs_joinids()))
  # Test named
  expect_no_condition(named <- query$to_seurat_assay(
    cells_index = 'baz',
    features_index = 'quux'
  ))
  expect_identical(
    rownames(named),
    query$var('quux')$GetColumnByName('quux')$as_vector()
  )
  expect_identical(
    colnames(named),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
})

test_that("Load assay from indexed ExperimentQuery", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("soma-experiment-query-value-filters")
  n_obs <- 1001L
  n_var <- 99L
  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )
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
  expect_no_condition(assay <- query$to_seurat_assay())
  expect_s4_class(assay, 'Assay')
  expect_identical(
    dim(assay),
    c(length(var_label_values), length(obs_label_values))
  )
  expect_identical(rownames(assay), paste0('feature', query$var_joinids()))
  expect_identical(colnames(assay), paste0('cell', query$obs_joinids()))
  # Test named
  expect_no_condition(named <- query$to_seurat_assay(
    cells_index = 'baz',
    features_index = 'quux'
  ))
  expect_identical(
    rownames(named),
    query$var('quux')$GetColumnByName('quux')$as_vector()
  )
  expect_identical(rownames(named), var_label_values)
  expect_identical(
    colnames(named),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_identical(colnames(named), obs_label_values)
})

test_that("Load reduction from ExperimentQuery mechanics", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("reduc-experiment-query-whole")
  n_obs <- 20L
  n_var <- 10L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )
  # Add embeddings
  n_pcs <- 50L
  n_umaps <- 2L
  obsm <- SOMACollection$new(file.path(experiment$ms$get('RNA')$uri, 'obsm'))
  obsm$create()
  obsm$add_new_dense_ndarray(
    key = 'X_pca',
    type = arrow::int32(),
    shape = c(n_obs, n_pcs)
  )
  obsm$get('X_pca')$write(create_dense_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_pcs
  ))
  obsm$add_new_dense_ndarray(
    key = 'X_umap',
    type = arrow::int32(),
    shape = c(n_obs, n_umaps)
  )
  obsm$get('X_umap')$write(create_dense_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_umaps,
    seed = 2L
  ))
  experiment$ms$get("RNA")$add_new_collection(obsm, 'obsm')
  # Add loadings
  varm <- SOMACollection$new(file.path(experiment$ms$get('RNA')$uri, 'varm'))
  varm$create()
  varm$add_new_dense_ndarray(
    key = 'PCs',
    type = arrow::int32(),
    shape = c(n_var, n_pcs)
  )
  varm$get('PCs')$write(create_dense_matrix_with_int_dims(
    nrows = n_var,
    ncols = n_pcs
  ))
  experiment$ms$get('RNA')$add_new_collection(varm, 'varm')
  # Create the query
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )
  # Test loading reductions
  expect_no_condition(X_pca <- query$to_seurat_reduction('X_pca'))
  expect_s4_class(X_pca, 'DimReduc')
  expect_equal(dim(X_pca), c(n_obs, n_pcs))
  expect_equal(dim(SeuratObject::Loadings(X_pca)), c(n_var, n_pcs))
  expect_false(SeuratObject::IsGlobal(X_pca))
  expect_equal(SeuratObject::Key(X_pca), 'PC_')
  expect_no_condition(X_umap <- query$to_seurat_reduction('X_umap'))
  expect_s4_class(X_umap, 'DimReduc')
  expect_equal(dim(X_umap), c(n_obs, n_umaps))
  expect_equal(dim(SeuratObject::Loadings(X_umap)), c(0L, 0L))
  expect_true(SeuratObject::IsGlobal(X_umap))
  expect_equal(SeuratObject::Key(X_umap), 'UMAP_')
  # Test using Seurat names
  expect_no_condition(pca <- query$to_seurat_reduction('pca'))
  expect_s4_class(pca, 'DimReduc')
  expect_true(is.matrix(SeuratObject::Embeddings(pca)))
  expect_true(is.matrix(SeuratObject::Loadings(pca)))
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(pca, TRUE)))
  expect_false(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(pca, FALSE)))
  expect_equal(dim(pca), c(n_obs, n_pcs))
  expect_equal(dim(SeuratObject::Embeddings(pca)), dim(pca))
  expect_equal(dim(SeuratObject::Loadings(pca)), c(n_var, n_pcs))
  expect_equal(dim(SeuratObject::Loadings(pca, TRUE)), c(0L, 0L))
  expect_false(SeuratObject::IsGlobal(pca))
  expect_equal(SeuratObject::Key(pca), 'PC_')
  expect_identical(
    colnames(SeuratObject::Embeddings(pca)),
    paste0('PC_', seq_len(n_pcs))
  )
  expect_identical(
    colnames(SeuratObject::Loadings(pca)),
    paste0('PC_', seq_len(n_pcs))
  )
  expect_no_condition(umap <- query$to_seurat_reduction('umap'))
  expect_s4_class(umap, 'DimReduc')
  expect_true(is.matrix(SeuratObject::Embeddings(umap)))
  expect_true(is.matrix(SeuratObject::Loadings(umap)))
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(umap, TRUE)))
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(umap, FALSE)))
  expect_equal(dim(umap), c(n_obs, n_umaps))
  expect_equal(dim(SeuratObject::Embeddings(umap)), dim(umap))
  expect_equal(dim(SeuratObject::Loadings(umap)), c(0L, 0L))
  expect_equal(dim(SeuratObject::Loadings(umap, TRUE)), c(0L, 0L))
  expect_true(SeuratObject::IsGlobal(umap))
  expect_equal(SeuratObject::Key(umap), 'UMAP_')
  expect_identical(
    colnames(SeuratObject::Embeddings(umap)),
    paste0('UMAP_', seq_len(n_umaps))
  )
  # Test adding names
  expect_no_condition(named_pca <- query$to_seurat_reduction(
    'pca',
    cells_index = 'baz',
    features_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_pca),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Embeddings(named_pca)),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Loadings(named_pca)),
    query$var('quux')$GetColumnByName('quux')$as_vector()
  )
  expect_no_condition(named_umap <- query$to_seurat_reduction(
    'umap',
    cells_index = 'baz',
  ))
  expect_identical(
    SeuratObject::Cells(named_umap),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Embeddings(named_umap)),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_no_condition(named_umap <- query$to_seurat_reduction(
    'umap',
    cells_index = 'baz',
    features_index = 'quux'
  ))
  # Test suppressing feature loadings
  suppress <- list(NA, FALSE)
  for (i in seq_along(suppress)) {
    expect_no_condition(no_load <- query$to_seurat_reduction(
      'pca',
      suppress[[i]]
    ))
    expect_equal(dim(no_load), c(n_obs, n_pcs))
    expect_equal(dim(SeuratObject::Loadings(no_load)), c(0L, 0L))
    expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(no_load)))
  }
  # Test assertions
  expect_error(query$to_seurat_reduction(TRUE))
  expect_error(query$to_seurat_reduction(NULL))
  expect_error(query$to_seurat_reduction(1))
  expect_error(query$to_seurat_reduction(c('pca', 'umap')))
  expect_error(query$to_seurat_reduction('tomato'))
  expect_error(query$to_seurat_reduction('pca', 1))
  expect_error(query$to_seurat_reduction('pca', cells_index = FALSE))
  expect_error(query$to_seurat_reduction('pca', cells_index = NA_character_))
  expect_error(query$to_seurat_reduction('pca', cells_index = 1))
  expect_error(query$to_seurat_reduction('pca', cells_index = c('baz', 'foo')))
  expect_error(query$to_seurat_reduction('pca', cells_index = 'tomato'))
  # Test `features_index` assertions
  expect_error(query$to_seurat_reduction('pca', features_index = FALSE))
  expect_error(query$to_seurat_reduction('pca', features_index = NA_character_))
  expect_error(query$to_seurat_reduction('pca', features_index = 1))
  expect_error(query$to_seurat_reduction(
    'pca',
    features_index = c('baz', 'foo')
  ))
  expect_error(query$to_seurat_reduction('pca', features_index = 'tomato'))
})

test_that("Load reduction from sliced ExperimentQuery", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("reduction-experiment-query-sliced")
  n_obs <- 1001L
  n_var <- 99L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )
  # Add embeddings
  n_pcs <- 50L
  n_umaps <- 2L
  obsm <- SOMACollection$new(file.path(experiment$ms$get('RNA')$uri, 'obsm'))
  obsm$create()
  obsm$add_new_dense_ndarray(
    key = 'X_pca',
    type = arrow::int32(),
    shape = c(n_obs, n_pcs)
  )
  obsm$get('X_pca')$write(create_dense_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_pcs
  ))
  obsm$add_new_dense_ndarray(
    key = 'X_umap',
    type = arrow::int32(),
    shape = c(n_obs, n_umaps)
  )
  obsm$get('X_umap')$write(create_dense_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_umaps,
    seed = 2L
  ))
  experiment$ms$get("RNA")$add_new_collection(obsm, 'obsm')
  # Add loadings
  varm <- SOMACollection$new(file.path(experiment$ms$get('RNA')$uri, 'varm'))
  varm$create()
  varm$add_new_dense_ndarray(
    key = 'PCs',
    type = arrow::int32(),
    shape = c(n_var, n_pcs)
  )
  varm$get('PCs')$write(create_dense_matrix_with_int_dims(
    nrows = n_var,
    ncols = n_pcs
  ))
  experiment$ms$get('RNA')$add_new_collection(varm, 'varm')
  # Create the query
  obs_slice <- bit64::as.integer64(seq(3, 72))
  var_slice <- bit64::as.integer64(seq(7, 21))
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(coords = list(soma_joinid = obs_slice)),
    var_query = SOMAAxisQuery$new(coords = list(soma_joinid = var_slice))
  )
  expect_no_condition(pca <- query$to_seurat_reduction('pca'))
  expect_s4_class(pca, 'DimReduc')
  expect_identical(dim(pca), c(length(obs_slice), n_pcs))
  expect_identical(dim(SeuratObject::Embeddings(pca)), dim(pca))
  expect_identical(dim(SeuratObject::Loadings(pca)), c(length(var_slice), n_pcs))
  expect_identical(SeuratObject::Cells(pca), paste0('cell', query$obs_joinids()))
  expect_identical(
    rownames(SeuratObject::Loadings(pca)),
    paste0('feature', query$var_joinids())
  )
  expect_false(SeuratObject::IsGlobal(pca))
  expect_equal(SeuratObject::Key(pca), 'PC_')
  expect_identical(
    colnames(SeuratObject::Embeddings(pca)),
    paste0('PC_', seq_len(n_pcs))
  )
  expect_identical(
    colnames(SeuratObject::Loadings(pca)),
    paste0('PC_', seq_len(n_pcs))
  )
  expect_no_condition(umap <- query$to_seurat_reduction('umap'))
  expect_s4_class(umap, 'DimReduc')
  expect_identical(dim(umap), c(length(obs_slice), n_umaps))
  expect_identical(dim(SeuratObject::Embeddings(umap)), dim(umap))
  expect_identical(dim(SeuratObject::Loadings(umap)), c(0L, 0L))
  expect_identical(
    SeuratObject::Cells(umap),
    paste0('cell', query$obs_joinids())
  )
  expect_true(SeuratObject::IsGlobal(umap))
  expect_equal(SeuratObject::Key(umap), 'UMAP_')
  expect_identical(
    colnames(SeuratObject::Embeddings(umap)),
    paste0('UMAP_', seq_len(n_umaps))
  )
  # Test named
  expect_no_condition(named_pca <- query$to_seurat_reduction(
    'pca',
    cells_index = 'baz',
    features_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_pca),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Loadings(named_pca)),
    query$var('quux')$GetColumnByName('quux')$as_vector()
  )
  expect_no_condition(named_umap <- query$to_seurat_reduction(
    'umap',
    cells_index = 'baz'
  ))
  expect_identical(
    SeuratObject::Cells(named_umap),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_no_condition(query$to_seurat_reduction(
    'umap',
    cells_index = 'baz',
    features_index = 'quux'
  ))
})

test_that("Load reduction from indexed ExperimentQuery", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("reduction-experiment-query-value-filters")
  n_obs <- 1001L
  n_var <- 99L
  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts")
  )
  # Add embeddings
  n_pcs <- 50L
  n_umaps <- 2L
  obsm <- SOMACollection$new(file.path(experiment$ms$get('RNA')$uri, 'obsm'))
  obsm$create()
  obsm$add_new_dense_ndarray(
    key = 'X_pca',
    type = arrow::int32(),
    shape = c(n_obs, n_pcs)
  )
  obsm$get('X_pca')$write(create_dense_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_pcs
  ))
  obsm$add_new_dense_ndarray(
    key = 'X_umap',
    type = arrow::int32(),
    shape = c(n_obs, n_umaps)
  )
  obsm$get('X_umap')$write(create_dense_matrix_with_int_dims(
    nrows = n_obs,
    ncols = n_umaps,
    seed = 2L
  ))
  experiment$ms$get("RNA")$add_new_collection(obsm, 'obsm')
  # Add loadings
  varm <- SOMACollection$new(file.path(experiment$ms$get('RNA')$uri, 'varm'))
  varm$create()
  varm$add_new_dense_ndarray(
    key = 'PCs',
    type = arrow::int32(),
    shape = c(n_var, n_pcs)
  )
  varm$get('PCs')$write(create_dense_matrix_with_int_dims(
    nrows = n_var,
    ncols = n_pcs
  ))
  experiment$ms$get('RNA')$add_new_collection(varm, 'varm')
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
  expect_no_condition(pca <- query$to_seurat_reduction('pca'))
  expect_s4_class(pca, 'DimReduc')
  expect_identical(
    dim(pca),
    c(length(obs_label_values), n_pcs)
  )
  expect_identical(dim(SeuratObject::Embeddings(pca)), dim(pca))
  expect_identical(
    dim(SeuratObject::Loadings(pca)),
    c(length(var_label_values), n_pcs)
  )
  expect_identical(SeuratObject::Cells(pca), paste0('cell', query$obs_joinids()))
  expect_identical(
    rownames(SeuratObject::Loadings(pca)),
    paste0('feature', query$var_joinids())
  )
  expect_false(SeuratObject::IsGlobal(pca))
  expect_equal(SeuratObject::Key(pca), 'PC_')
  expect_identical(
    colnames(SeuratObject::Embeddings(pca)),
    paste0('PC_', seq_len(n_pcs))
  )
  expect_identical(
    colnames(SeuratObject::Loadings(pca)),
    paste0('PC_', seq_len(n_pcs))
  )
  expect_no_condition(umap <- query$to_seurat_reduction('umap'))
  expect_s4_class(umap, 'DimReduc')
  expect_identical(
    dim(umap),
    c(length(obs_label_values), n_umaps)
  )
  expect_identical(dim(SeuratObject::Embeddings(umap)), dim(umap))
  expect_identical(dim(SeuratObject::Loadings(umap)), c(0L, 0L))
  expect_identical(
    SeuratObject::Cells(umap),
    paste0('cell', query$obs_joinids())
  )
  expect_true(SeuratObject::IsGlobal(umap))
  expect_equal(SeuratObject::Key(umap), 'UMAP_')
  expect_identical(
    colnames(SeuratObject::Embeddings(umap)),
    paste0('UMAP_', seq_len(n_umaps))
  )
  # Test named
  expect_no_condition(named_pca <- query$to_seurat_reduction(
    'pca',
    cells_index = 'baz',
    features_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_pca),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::Cells(named_pca), obs_label_values)
  expect_identical(
    rownames(SeuratObject::Loadings(named_pca)),
    query$var('quux')$GetColumnByName('quux')$as_vector()
  )
  expect_identical(rownames(SeuratObject::Loadings(named_pca)), var_label_values)
  expect_no_condition(named_umap <- query$to_seurat_reduction(
    'umap',
    cells_index = 'baz'
  ))
  expect_identical(
    SeuratObject::Cells(named_umap),
    query$obs('baz')$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::Cells(named_umap), obs_label_values)
  expect_no_condition(query$to_seurat_reduction(
    'umap',
    cells_index = 'baz',
    features_index = 'quux'
  ))
})
