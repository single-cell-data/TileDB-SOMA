test_that("Load reduction from ExperimentQuery mechanics", {
  skip_if(!extended_tests()  || covr_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("reduc-experiment-query-whole")
  n_obs <- 20L
  n_var <- 10L
  n_pcs <- 50L
  n_ics <- 30L
  n_umaps <- 2L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    obsm_layers = c(X_pca = n_pcs, "dense:X_ica" = n_ics, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs, "dense:ICs" = n_ics),
    mode = "READ"
  )
  on.exit(experiment$close())

  # Create the query
  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  # Test loading reductions
  expect_no_condition(X_pca <- query$to_seurat_reduction('X_pca'))
  expect_s4_class(X_pca, 'DimReduc')
  expect_equal(dim(X_pca), c(n_obs, n_pcs))

  expect <- dim(SeuratObject::Loadings(X_pca))
  actual <- c(n_var, n_pcs)

  expect_equal(expect, actual)

  expect_false(SeuratObject::IsGlobal(X_pca))
  expect_equal(SeuratObject::Key(X_pca), 'PC_')
  expect_warning(X_ica <- query$to_seurat_reduction('X_ica'))
  expect_s4_class(X_ica, 'DimReduc')
  expect_equal(dim(X_ica), c(n_obs, n_ics))
  expect_false(SeuratObject::IsGlobal(X_ica))
  expect_equal(SeuratObject::Key(X_ica), 'IC_')
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
  expect_warning(ica <- query$to_seurat_reduction('ica'))
  expect_s4_class(ica, 'DimReduc')
  expect_true(is.matrix(SeuratObject::Embeddings(ica)))
  expect_true(is.matrix(SeuratObject::Loadings(ica)))
  expect_true(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(ica, TRUE)))
  expect_false(SeuratObject::IsMatrixEmpty(SeuratObject::Loadings(ica, FALSE)))
  expect_equal(dim(ica), c(n_obs, n_ics))
  expect_equal(dim(SeuratObject::Embeddings(ica)), dim(ica))
  expect_equal(dim(SeuratObject::Loadings(ica)), c(n_var, n_ics))
  expect_equal(dim(SeuratObject::Loadings(ica, TRUE)), c(0L, 0L))
  expect_false(SeuratObject::IsGlobal(ica))
  expect_equal(SeuratObject::Key(ica), 'IC_')
  expect_identical(
    colnames(SeuratObject::Embeddings(ica)),
    paste0('IC_', seq_len(n_ics))
  )
  expect_identical(
    colnames(SeuratObject::Loadings(ica)),
    paste0('IC_', seq_len(n_ics))
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
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_pca),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Embeddings(named_pca)),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Loadings(named_pca)),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_warning(named_ica <- query$to_seurat_reduction(
    'ica',
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_ica),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Embeddings(named_ica)),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Loadings(named_ica)),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_no_condition(named_umap <- query$to_seurat_reduction(
    'umap',
    obs_index = 'baz',
  ))
  expect_identical(
    SeuratObject::Cells(named_umap),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Embeddings(named_umap)),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_no_condition(named_umap <- query$to_seurat_reduction(
    'umap',
    obs_index = 'baz',
    var_index = 'quux'
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
  expect_error(query$to_seurat_reduction('pca', 'LOADINGS'))
  expect_error(query$to_seurat_reduction('pca', obs_index = FALSE))
  expect_error(query$to_seurat_reduction('pca', obs_index = NA_character_))
  expect_error(query$to_seurat_reduction('pca', obs_index = 1))
  expect_error(query$to_seurat_reduction('pca', obs_index = c('baz', 'foo')))
  expect_error(query$to_seurat_reduction('pca', obs_index = 'tomato'))

  # Test `var_index` assertions
  expect_error(query$to_seurat_reduction('pca', var_index = FALSE))
  expect_error(query$to_seurat_reduction('pca', var_index = NA_character_))
  expect_error(query$to_seurat_reduction('pca', var_index = 1))
  expect_error(query$to_seurat_reduction(
    'pca',
    var_index = c('baz', 'foo')
  ))
  expect_error(query$to_seurat_reduction('pca', var_index = 'tomato'))
})

test_that("Load reduction from sliced ExperimentQuery", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("reduction-experiment-query-sliced")
  n_obs <- 1001L
  n_var <- 99L
  n_pcs <- 50L
  n_ics <- 30L
  n_umaps <- 2L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    obsm_layers = c(X_pca = n_pcs, "dense:X_ica" = n_ics, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs, "dense:ICs" = n_ics),
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
  expect_no_condition(pca <- query$to_seurat_reduction('pca'))
  expect_s4_class(pca, 'DimReduc')
  expect_identical(dim(pca), c(length(obs_slice), n_pcs))
  expect_identical(dim(SeuratObject::Embeddings(pca)), dim(pca))
  expect_identical(dim(SeuratObject::Loadings(pca)), c(length(var_slice), n_pcs))
  expect_identical(SeuratObject::Cells(pca), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(
    rownames(SeuratObject::Loadings(pca)),
    paste0('feature', query$var_joinids()$as_vector())
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
  expect_warning(ica <- query$to_seurat_reduction('ica'))
  expect_s4_class(ica, 'DimReduc')
  expect_identical(dim(ica), c(length(obs_slice), n_ics))
  expect_identical(dim(SeuratObject::Embeddings(ica)), dim(ica))
  expect_identical(dim(SeuratObject::Loadings(ica)), c(length(var_slice), n_ics))
  expect_identical(SeuratObject::Cells(ica), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(
    rownames(SeuratObject::Loadings(ica)),
    paste0('feature', query$var_joinids()$as_vector())
  )
  expect_false(SeuratObject::IsGlobal(ica))
  expect_equal(SeuratObject::Key(ica), 'IC_')
  expect_identical(
    colnames(SeuratObject::Embeddings(ica)),
    paste0('IC_', seq_len(n_ics))
  )
  expect_identical(
    colnames(SeuratObject::Loadings(ica)),
    paste0('IC_', seq_len(n_ics))
  )
  expect_no_condition(umap <- query$to_seurat_reduction('umap'))
  expect_s4_class(umap, 'DimReduc')
  expect_identical(dim(umap), c(length(obs_slice), n_umaps))
  expect_identical(dim(SeuratObject::Embeddings(umap)), dim(umap))
  expect_identical(dim(SeuratObject::Loadings(umap)), c(0L, 0L))
  expect_identical(
    SeuratObject::Cells(umap),
    paste0('cell', query$obs_joinids()$as_vector())
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
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_pca),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Loadings(named_pca)),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_warning(named_ica <- query$to_seurat_reduction(
    'ica',
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_ica),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(SeuratObject::Loadings(named_ica)),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_no_condition(named_umap <- query$to_seurat_reduction(
    'umap',
    obs_index = 'baz'
  ))
  expect_identical(
    SeuratObject::Cells(named_umap),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_no_condition(query$to_seurat_reduction(
    'umap',
    obs_index = 'baz',
    var_index = 'quux'
  ))
})

test_that("Load reduction from indexed ExperimentQuery", {
  skip_if(!extended_tests() || covr_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- withr::local_tempdir("reduction-experiment-query-value-filters")
  n_obs <- 1001L
  n_var <- 99L
  n_pcs <- 50L
  n_ics <- 30L
  n_umaps <- 2L
  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    obsm_layers = c(X_pca = n_pcs, "dense:X_ica" = n_ics, X_umap = n_umaps),
    varm_layers = c(PCs = n_pcs, "dense:ICs" = n_ics),
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
  expect_identical(SeuratObject::Cells(pca), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(
    rownames(SeuratObject::Loadings(pca)),
    paste0('feature', query$var_joinids()$as_vector())
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
  expect_warning(ica <- query$to_seurat_reduction('ica'))
  expect_s4_class(ica, 'DimReduc')
  expect_identical(
    dim(ica),
    c(length(obs_label_values), n_ics)
  )
  expect_identical(dim(SeuratObject::Embeddings(ica)), dim(ica))
  expect_identical(
    dim(SeuratObject::Loadings(ica)),
    c(length(var_label_values), n_ics)
  )
  expect_identical(SeuratObject::Cells(ica), paste0('cell', query$obs_joinids()$as_vector()))
  expect_identical(
    rownames(SeuratObject::Loadings(ica)),
    paste0('feature', query$var_joinids()$as_vector())
  )
  expect_false(SeuratObject::IsGlobal(ica))
  expect_equal(SeuratObject::Key(ica), 'IC_')
  expect_identical(
    colnames(SeuratObject::Embeddings(ica)),
    paste0('IC_', seq_len(n_ics))
  )
  expect_identical(
    colnames(SeuratObject::Loadings(ica)),
    paste0('IC_', seq_len(n_ics))
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
    paste0('cell', query$obs_joinids()$as_vector())
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
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_pca),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::Cells(named_pca), obs_label_values)
  expect_identical(
    rownames(SeuratObject::Loadings(named_pca)),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_identical(rownames(SeuratObject::Loadings(named_pca)), var_label_values)
  expect_warning(named_ica <- query$to_seurat_reduction(
    'ica',
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    SeuratObject::Cells(named_ica),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::Cells(named_ica), obs_label_values)
  expect_identical(
    rownames(SeuratObject::Loadings(named_ica)),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_identical(rownames(SeuratObject::Loadings(named_ica)), var_label_values)
  expect_no_condition(named_umap <- query$to_seurat_reduction(
    'umap',
    obs_index = 'baz'
  ))
  expect_identical(
    SeuratObject::Cells(named_umap),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(SeuratObject::Cells(named_umap), obs_label_values)
  expect_no_condition(query$to_seurat_reduction(
    'umap',
    obs_index = 'baz',
    var_index = 'quux'
  ))
})
