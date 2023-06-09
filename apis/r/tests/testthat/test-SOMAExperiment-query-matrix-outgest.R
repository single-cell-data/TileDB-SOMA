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
  # Need to manually name the columns since we don't store them
  colnames(pcs2) <- paste0("PC_", 1:ncol(pcs2))
  expect_identical(pcs1, as.matrix(pcs2[rownames(pcs1), colnames(pcs2)]))

  # Embeddings
  pcas1 <- SeuratObject::Embeddings(pbmc_small, "pca")
  pcas2 <- query$to_matrix(
    collection = "obsm",
    layer_name = "X_pca",
    obs_index = "obs_id"
  )
  colnames(pcas2) <- paste0("PC_", 1:ncol(pcas2))
  expect_identical(pcas1, as.matrix(pcas2[rownames(pcas1), colnames(pcas1)]))

  # Graphs
  # Need to coerce original and retrieved graphs to dgCMatrices for comparison
  snn1 <- as(SeuratObject::Graphs(pbmc_small, "RNA_snn"), "CsparseMatrix")
  snn2 <- query$to_matrix(
    collection = "obsp",
    layer_name = "RNA_snn",
    obs_index = "obs_id"
  ) |> as("CsparseMatrix")
  expect_identical(snn1, snn2)


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
