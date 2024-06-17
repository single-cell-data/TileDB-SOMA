test_that("Load assay from ExperimentQuery mechanics", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))

  uri <- tempfile(pattern="assay-experiment-query-whole")
  n_obs <- 20L
  n_var <- 10L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

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
  expect_equal(rownames(assay), paste0('feature', query$var_joinids()$as_vector()))
  expect_equal(colnames(assay), paste0('cell', seq_len(n_obs) - 1L))
  expect_equal(colnames(assay), paste0('cell', query$obs_joinids()$as_vector()))

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
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    colnames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(
    rownames(named),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
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
  expect_error(query$to_seurat_assay(c(counts = 'counts', data = 'tomato')))

  # Test `obs_index` assertions
  expect_error(query$to_seurat_assay(obs_index = FALSE))
  expect_error(query$to_seurat_assay(obs_index = NA_character_))
  expect_error(query$to_seurat_assay(obs_index = 1))
  expect_error(query$to_seurat_assay(obs_index = c('baz', 'foo')))
  expect_error(query$to_seurat_assay(obs_index = 'tomato'))

  # Test `var_index` assertions
  expect_error(query$to_seurat_assay(var_index = FALSE))
  expect_error(query$to_seurat_assay(var_index = NA_character_))
  expect_error(query$to_seurat_assay(var_index = 1))
  expect_error(query$to_seurat_assay(var_index = c('baz', 'foo')))
  expect_error(query$to_seurat_assay(var_index = 'tomato'))

  # Test `var_column_names` assertions
  expect_error(query$to_seurat_assay(var_column_names = 1L))
  expect_error(query$to_seurat_assay(var_column_names = c(
    NA_character_,
    NA_character_
  )))
  expect_error(query$to_seurat_assay(var_column_names = c(TRUE, FALSE)))
  expect_error(query$to_seurat_assay(var_column_names = 'tomato'))
})

test_that("Load assay with SeuratObject v5 returns v3 assays", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', '4.9.9.9094')

  withr::local_options(Seurat.object.assay.version = 'v5')
  uri <- tempfile(pattern="assay-experiment-query-v5-v3")
  n_obs <- 20L
  n_var <- 10L
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  expect_identical(getOption('Seurat.object.assay.version'), 'v5')
  expect_no_condition(assay <- query$to_seurat_assay())
  expect_identical(getOption('Seurat.object.assay.version'), 'v5')
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
  expect_equal(rownames(assay), paste0('feature', query$var_joinids()$as_vector()))
  expect_equal(colnames(assay), paste0('cell', seq_len(n_obs) - 1L))
  expect_equal(colnames(assay), paste0('cell', query$obs_joinids()$as_vector()))
})

test_that("Load assay from sliced ExperimentQuery", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))

  uri <- tempfile(pattern="assay-experiment-query-sliced")
  n_obs <- 1001L
  n_var <- 99L
  obs_slice <- bit64::as.integer64(seq(3, 72))
  var_slice <- bit64::as.integer64(seq(7, 21))
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    obs_query = SOMAAxisQuery$new(coords = list(soma_joinid = obs_slice)),
    var_query = SOMAAxisQuery$new(coords = list(soma_joinid = var_slice))
  )
  expect_no_condition(assay <- query$to_seurat_assay())
  expect_s4_class(assay, 'Assay')
  expect_identical(dim(assay), c(length(var_slice), length(obs_slice)))
  expect_identical(rownames(assay), paste0('feature', query$var_joinids()$as_vector()))
  expect_identical(colnames(assay), paste0('cell', query$obs_joinids()$as_vector()))

  # Test named
  expect_no_condition(named <- query$to_seurat_assay(
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    rownames(named),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_identical(
    colnames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
})

test_that("Load assay from indexed ExperimentQuery", {
  skip_if(!extended_tests())
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  uri <- tempfile(pattern="soma-experiment-query-value-filters")
  n_obs <- 1001L
  n_var <- 99L
  obs_label_values <- c("1003", "1007", "1038", "1099")
  var_label_values <- c("1018", "1034", "1067")
  experiment <- create_and_populate_experiment(
    uri = uri,
    n_obs = n_obs,
    n_var = n_var,
    X_layer_names = c("counts", "logcounts"),
    mode = "READ"
  )
  on.exit(experiment$close())

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
  expect_identical(rownames(assay), paste0('feature', query$var_joinids()$as_vector()))
  expect_identical(colnames(assay), paste0('cell', query$obs_joinids()$as_vector()))

  # Test named
  expect_no_condition(named <- query$to_seurat_assay(
    obs_index = 'baz',
    var_index = 'quux'
  ))
  expect_identical(
    rownames(named),
    query$var('quux')$concat()$GetColumnByName('quux')$as_vector()
  )
  expect_identical(rownames(named), var_label_values)
  expect_identical(
    colnames(named),
    query$obs('baz')$concat()$GetColumnByName('baz')$as_vector()
  )
  expect_identical(colnames(named), obs_label_values)
})
