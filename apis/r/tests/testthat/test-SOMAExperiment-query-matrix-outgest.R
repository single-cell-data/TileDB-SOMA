test_that("matrix outgest with all results", {
  experiment <- load_dataset("soma-exp-pbmc-small")

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  # Loadings
  pcs1 <- SeuratObject::Loadings(pbmc_small, "pca")
  pcs2 <- query$to_matrix(
    collection = "varm",
    layer_name = "PCs",
    var_index = "var_id"
  )
  # Column names are unset since we don't store them
  expect_null(colnames(pcs2))
  # Manually name the columns for comparison
  colnames(pcs2) <- paste0("PC_", 1:ncol(pcs2))
  # Coerce to dense matrix and reorder for comparison
  as.matrix(pcs2[rownames(pcs1), colnames(pcs2)])
  expect_identical(pcs2, pcs1)

  # Embeddings
  pcas1 <- SeuratObject::Embeddings(pbmc_small, "pca")
  pcas2 <- query$to_matrix(
    collection = "obsm",
    layer_name = "X_pca",
    obs_index = "obs_id"
  )
  expect_null(colnames(pcs2))
  colnames(pcas2) <- paste0("PC_", 1:ncol(pcas2))
  pcas2 <- as.matrix(pcas2[rownames(pcas1), colnames(pcas1)])
  expect_identical(pcas2, pcas1)

  # Graphs
  # Need to coerce original and retrieved graphs to dgCMatrices for comparison
  snn1 <- as(SeuratObject::Graphs(pbmc_small, "RNA_snn"), "CsparseMatrix")
  snn2 <- query$to_matrix(
    collection = "obsp",
    layer_name = "RNA_snn",
    obs_index = "obs_id"
  ) |> as("CsparseMatrix")
  expect_identical(snn1, snn2)

  # Assay data
  assay1 <- SeuratObject::GetAssayData(pbmc_small[["RNA"]], slot = "counts")
  assay2 <- query$to_matrix(
    collection = "X",
    layer_name = "counts",
    obs_index = "obs_id",
    var_index = "var_id"
  )

  # Transpose and coerce to dGCMatrix for comparison
  assay2 <- as(t(assay2), "CsparseMatrix")
  expect_equal(assay2, assay1)
})

test_that("matrix outgest assertions", {
  experiment <- load_dataset("soma-exp-pbmc-small")

  query <- SOMAExperimentAxisQuery$new(
    experiment = experiment,
    measurement_name = "RNA"
  )

  expect_error(
    query$to_matrix(collection = "foo"),
    "The following collection does not exist: foo"
  )

  expect_error(
    query$to_matrix(collection = "X", layer_name = "foo"),
    "The following layer does not exist: foo"
  )

  expect_message(
    query$to_matrix(collection = "X", layer_name = "counts", obs_index = "foo"),
    "The following obs index does not exist: foo"
  )

  # an unnamed matrix is returned if obs/var_index args are not provided
  expect_identical(
    dimnames(query$to_matrix("X", "counts")),
    list(NULL, NULL)
  )

  # error if specified obs/var_index does not exist
  expect_error(
    query$to_matrix("X", "counts", obs_index = "foo"),
    "The following column does not exist: foo"
  )
  expect_error(
    query$to_matrix("X", "counts", var_index = "foo"),
    "The following column does not exist: foo"
  )

  # only one of obs_index or var_index can be provided
  expect_identical(
    dimnames(query$to_matrix("X", "counts", obs_index = "obs_id")),
    list(query$obs(column_names = "obs_id")$obs_id$as_vector(), NULL)
  )
  expect_identical(
    dimnames(query$to_matrix("X", "counts", var_index = "var_id")),
    list(NULL, query$var(column_names = "var_id")$var_id$as_vector())
  )

  # a warning is issued if an index is unnecessarily provided
  expect_output(
    query$to_matrix("obsm", "X_pca", var_index = "var_id"),
    "The var_index is ignored for obsm collections"
  )
  expect_output(
    query$to_matrix("varm", "PCs", obs_index = "obs_id"),
    "The obs_index is ignored for varm collections"
  )

  # retrieved labels must be unique
  expect_error(
    query$to_matrix("obsm", "X_pca", obs_index = "groups"),
    "All obs_index values must be unique"
  )
})
