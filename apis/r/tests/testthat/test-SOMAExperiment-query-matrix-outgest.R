test_that("matrix outgest with all results", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  pbmc_small <- get_data("pbmc_small", package = "SeuratObject")
  experiment <- load_dataset("soma-exp-pbmc-small")

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  # Loadings
  pcs1 <- SeuratObject::Loadings(pbmc_small, "pca")
  pcs2 <- query$to_sparse_matrix(
    collection = "varm",
    layer_name = "PCs",
    var_index = "var_id"
  )

  # Column names for non obs/var dimensions default to soma_dim_1 values
  expect_equal(colnames(pcs2), as.character(seq_len(ncol(pcs2)) - 1))

  # Manually name the columns for comparison
  colnames(pcs2) <- paste0("PC_", 1:ncol(pcs2))
  # Coerce to dense matrix and reorder for comparison
  pcs2 <- as.matrix(pcs2[rownames(pcs1), colnames(pcs2)])
  expect_identical(pcs2, pcs1)

  # Embeddings
  pcas1 <- SeuratObject::Embeddings(pbmc_small, "pca")
  pcas2 <- query$to_sparse_matrix(
    collection = "obsm",
    layer_name = "X_pca",
    obs_index = "obs_id"
  )
  expect_equal(colnames(pcas2), as.character(seq_len(ncol(pcas2)) - 1))
  colnames(pcas2) <- paste0("PC_", 1:ncol(pcas2))
  pcas2 <- as.matrix(pcas2[rownames(pcas1), colnames(pcas1)])
  expect_identical(pcas2, pcas1)

  # Graphs
  # Need to coerce original and retrieved graphs to dgCMatrices for comparison
  snn1 <- as(SeuratObject::Graphs(pbmc_small, "RNA_snn"), "CsparseMatrix")
  snn2 <- query$to_sparse_matrix(
    collection = "obsp",
    layer_name = "RNA_snn",
    obs_index = "obs_id"
  ) |> as("CsparseMatrix")
  expect_identical(snn1, snn2)

  # Assay data
  assay1 <- SeuratObject::GetAssayData(pbmc_small[["RNA"]], slot = "counts")
  assay2 <- query$to_sparse_matrix(
    collection = "X",
    layer_name = "counts",
    obs_index = "obs_id",
    var_index = "var_id"
  )

  # Transpose and coerce to dGCMatrix for comparison
  assay2 <- as(Matrix::t(assay2), "CsparseMatrix")
  expect_equal(assay2, assay1)
})

test_that("matrix outgest with filtered results", {
  skip_if_not_installed('SeuratObject', .MINIMUM_SEURAT_VERSION('c'))
  # Subset the pbmc_small object to match the filtered results
  pbmc_small <- get_data("pbmc_small", package = "SeuratObject")
  pbmc_small1 <- pbmc_small[
    SeuratObject::VariableFeatures(pbmc_small),
    pbmc_small$nCount_RNA > 400
  ]

  experiment <- load_dataset("soma-exp-pbmc-small")

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA",
    var_query = SOMAAxisQuery$new(
      value_filter = "vst.variable == 'TRUE'"
    ),
    obs_query = SOMAAxisQuery$new(
      value_filter = "nCount_RNA > 400"
    )
  )

  expect_equal(query$n_obs, ncol(pbmc_small1))
  expect_equal(query$n_vars, nrow(pbmc_small1))

  # Loadings
  pcs1 <- SeuratObject::Loadings(pbmc_small1, "pca")
  pcs2 <- query$to_sparse_matrix(
    collection = "varm",
    layer_name = "PCs",
    var_index = "var_id"
  )

  # Manually name the columns for comparison
  colnames(pcs2) <- paste0("PC_", 1:ncol(pcs2))
  # Coerce to dense matrix and reorder for comparison
  pcs2 <- as.matrix(pcs2[rownames(pcs1), ])
  expect_identical(pcs2, pcs1)

  # Embeddings
  pcas1 <- SeuratObject::Embeddings(pbmc_small1, "pca")
  pcas2 <- query$to_sparse_matrix(
    collection = "obsm",
    layer_name = "X_pca",
    obs_index = "obs_id"
  )

  colnames(pcas2) <- paste0("PC_", 1:ncol(pcas2))
  pcas2 <- as.matrix(pcas2[rownames(pcas1), ])
  expect_identical(pcas2, pcas1)

  # Graphs
  # Seurat drops graphs from indexed objects so we need to subset manually
  snn1 <- as(SeuratObject::Graphs(pbmc_small, "RNA_snn"), "CsparseMatrix")
  snn1 <- snn1[colnames(pbmc_small1), colnames(pbmc_small1)]

  snn2 <- query$to_sparse_matrix(
    collection = "obsp",
    layer_name = "RNA_snn",
    obs_index = "obs_id"
  ) |> as("CsparseMatrix")
  expect_identical(snn2, snn1)

  # Assay data
  assay1 <- SeuratObject::GetAssayData(pbmc_small1[["RNA"]], slot = "counts")
  assay2 <- query$to_sparse_matrix(
    collection = "X",
    layer_name = "counts",
    obs_index = "obs_id",
    var_index = "var_id"
  )

  # Transpose and coerce to dGCMatrix for comparison
  assay2 <- as(Matrix::t(assay2), "CsparseMatrix")
  expect_identical(assay2[rownames(assay1), colnames(assay1)], assay1)
})

test_that("matrix outgest assertions", {
  experiment <- load_dataset("soma-exp-pbmc-small")

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  expect_error(
    query$to_sparse_matrix(collection = "foo"),
    "The following collection does not exist: foo"
  )

  expect_error(
    query$to_sparse_matrix(collection = "X", layer_name = "foo"),
    "The following layer does not exist: foo"
  )

  expect_error(
    query$to_sparse_matrix(collection = "X", layer_name = "counts", obs_index = "foo"),
    "The following column does not exist: foo"
  )

  # joinds are used as default dimnames if no obs/var_index is provided
  expect_identical(
    dimnames(query$to_sparse_matrix("X", "counts")),
    lapply(
      list(query$obs_joinids(), query$var_joinids()),
      function(x) as.character(x$as_vector())
    )
  )

  # error if specified obs/var_index does not exist
  expect_error(
    query$to_sparse_matrix("X", "counts", obs_index = "foo"),
    "The following column does not exist: foo"
  )
  expect_error(
    query$to_sparse_matrix("X", "counts", var_index = "foo"),
    "The following column does not exist: foo"
  )

  # only one of obs_index or var_index can be provided
  expect_identical(
    dimnames(query$to_sparse_matrix("X", "counts", obs_index = "obs_id")),
    list(
      as.character(query$obs(column_names = "obs_id")$concat()$obs_id),
      as.character(query$var_joinids()$as_vector())
    )
  )
  expect_identical(
    dimnames(query$to_sparse_matrix("X", "counts", var_index = "var_id")),
    list(
      as.character(query$obs_joinids()$as_vector()),
      query$var(column_names = "var_id")$concat()$var_id$as_vector()
    )
  )

  # a warning is issued if an index is unnecessarily provided
  expect_output(
    query$to_sparse_matrix("obsm", "X_pca", var_index = "var_id"),
    "The var_index is ignored for obsm collections"
  )
  expect_output(
    query$to_sparse_matrix("varm", "PCs", obs_index = "obs_id"),
    "The obs_index is ignored for varm collections"
  )

  # retrieved labels must be unique
  expect_error(
    query$to_sparse_matrix("obsm", "X_pca", obs_index = "groups"),
    "All obs_index values must be unique"
  )
})

test_that("matrix outgest with implicitly-stored axes", {
  uri <- withr::local_tempdir("matrix-implicit")
  set.seed(seed = 42L)
  n_obs <- 15L
  n_var <- 10L
  n_dims <- 5L
  m_key <- "X_dims"
  p_key <- "graph"
  shape <- rep_len(.Machine$integer.max - 1L, length.out = 2L)
  # Create our experiment
  experiment <- SOMAExperimentCreate(uri = uri)
  experiment$obs <- create_and_populate_obs(
    uri = file.path(experiment$uri, "obs"),
    nrows = n_obs
  )
  experiment$ms <- SOMACollectionCreate(file.path(experiment$uri, "ms"))
  ms_rna <- SOMAMeasurementCreate(file.path(experiment$ms$uri, "RNA"))
  experiment$ms$add_new_collection(object = ms_rna, key = "RNA")
  ms_rna$var <- create_and_populate_var(
    uri = file.path(ms_rna$uri, "var"),
    nrows = n_var
  )
  # Generate X matrix
  ms_rna$X <- SOMACollectionCreate(file.path(ms_rna$uri, "X"))
  X_counts <- Matrix::rsparsematrix(
    nrow = n_obs,
    ncol = n_var,
    density = 0.3,
    repr = "T"
  )
  X_counts[nrow(X_counts) %/% 2L, ] <- 0
  X_counts[, ncol(X_counts) %/% 2L] <- 0
  X_array <- SOMASparseNDArrayCreate(
    uri = file.path(ms_rna$X$uri, "counts"),
    type = arrow::infer_type(x = X_counts@x),
    shape = shape
  )
  X_array$write(X_counts)
  X_array$close()
  ms_rna$X$set(object = X_array, name = "counts")
  ms_rna$X$close()
  # Generate obsm
  ms_rna$obsm <- SOMACollectionCreate(file.path(ms_rna$uri, "obsm"))
  obsm <- Matrix::rsparsematrix(
    nrow = n_obs,
    ncol = n_dims,
    density = 0.7,
    repr = "T"
  )
  obsm[nrow(obsm) %/% 2L, ] <- 0
  obsm[, ncol(obsm) %/% 2L] <- 0
  ms_rna$obsm$add_new_sparse_ndarray(
    key = m_key,
    type = arrow::infer_type(x = obsm@x),
    shape = shape
  )
  ms_rna$obsm$get(m_key)$write(obsm)
  ms_rna$obsm$close()
  # Generate varm
  ms_rna$varm <- SOMACollectionCreate(file.path(ms_rna$uri, "varm"))
  varm <- Matrix::rsparsematrix(
    nrow = n_var,
    ncol = n_dims,
    density = 0.7,
    repr = "T"
  )
  varm[nrow(varm) %/% 2L, ] <- 0
  varm[, ncol(varm) %/% 2L] <- 0
  ms_rna$varm$add_new_sparse_ndarray(
    key = m_key,
    type = arrow::infer_type(x = varm@x),
    shape = shape
  )
  ms_rna$varm$get(m_key)$write(varm)
  ms_rna$varm$close()
  # Generate obsp
  ms_rna$obsp <- SOMACollectionCreate(file.path(ms_rna$uri, "obsp"))
  obsp <- Matrix::rsparsematrix(
    nrow = n_obs,
    ncol = n_obs,
    density = 0.3,
    repr = "T"
  )
  obsp[nrow(obsp) %/% 2L, ] <- 0
  obsp[, ncol(obsp) %/% 2L] <- 0
  ms_rna$obsp$add_new_sparse_ndarray(
    key = p_key,
    type = arrow::infer_type(x = obsp@x),
    shape = shape
  )
  ms_rna$obsp$get(p_key)$write(obsp)
  ms_rna$obsp$close()
  # Generate varp
  ms_rna$varp <- SOMACollectionCreate(file.path(ms_rna$uri, "varp"))
  varp <- Matrix::rsparsematrix(
    nrow = n_var,
    ncol = n_var,
    density = 0.3,
    repr = "T"
  )
  varp[nrow(varp) %/% 2L, ] <- 0
  varp[, ncol(varp) %/% 2L] <- 0
  ms_rna$varp$add_new_sparse_ndarray(
    key = p_key,
    type = arrow::infer_type(x = varp@x),
    shape = shape
  )
  ms_rna$varp$get(p_key)$write(varp)
  ms_rna$varp$close()
  # Close and reopen
  ms_rna$close()
  experiment$close()
  experiment <- SOMAExperimentOpen(experiment$uri)
  query <- SOMAExperimentAxisQuery$new(experiment, "RNA")
  # Read in matrices
  expect_s4_class(X_read <- query$to_sparse_matrix("X", "counts"), "dgTMatrix")
  expect_identical(dim(X_read), dim(X_counts))
  expect_s4_class(obsm_read <- query$to_sparse_matrix("obsm", m_key), "dgTMatrix")
  expect_identical(dim(obsm_read), dim(obsm))
  expect_s4_class(varm_read <- query$to_sparse_matrix("varm", m_key), "dgTMatrix")
  expect_identical(dim(varm_read), dim(varm))
  expect_s4_class(obsp_read <- query$to_sparse_matrix("obsp", p_key), "dgTMatrix")
  expect_identical(dim(obsp_read), dim(obsp))
  expect_s4_class(varp_read <- query$to_sparse_matrix("varp", p_key), "dgTMatrix")
  expect_identical(dim(varp_read), dim(varp))
})
