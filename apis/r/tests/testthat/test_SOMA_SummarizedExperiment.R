# TODO: Add tests for creating an SOMA from a SummarizedExperiment

test_that("a SummarizedExperiment can be created from an existing SOMA", {
  skip_if_not_installed("SummarizedExperiment")
  uri <- file.path(withr::local_tempdir(), "soma")

  assay <- pbmc_small[["RNA"]]
  soma <- SOMA$new(uri = uri)
  soma$from_seurat_assay(assay, obs = pbmc_small[[]])

  se_obj <- soma$to_summarized_experiment(layers = "counts")
  expect_s4_class(se_obj, "SummarizedExperiment")

  # use feature/sample names to ensure objects being compared are sorted
  var_ids <- rownames(assay)
  obs_ids <- colnames(assay)

  # validate sample metadata
  obs <- fac2char(pbmc_small[[]])

  expect_equal(
    as.data.frame(SummarizedExperiment::colData(se_obj))[obs_ids,],
    obs[obs_ids,]
  )

  # validate feature metadata
  # (manually remove vst.variable column because logicals are returned as ints)
  cols <- colnames(assay[[]])[-5]
  expect_equal(
    as.data.frame(SummarizedExperiment::rowData(se_obj))[var_ids, cols],
    assay[[]][var_ids, cols]
  )

  # validate raw counts matrix
  expect_identical(
    as.matrix(SummarizedExperiment::assays(se_obj)$counts[var_ids, obs_ids]),
    as.matrix(SeuratObject::GetAssayData(assay, "counts")[var_ids, obs_ids])
  )

  # validate normalized data matrix (and layer renaming)
  se_obj <- soma$to_summarized_experiment(layers = c(logcounts = "data"))
  expect_identical(
    as.matrix(SummarizedExperiment::assays(se_obj)$logcounts[var_ids, obs_ids]),
    as.matrix(SeuratObject::GetAssayData(assay, "data")[var_ids, obs_ids])
  )
})


test_that("a SingleCellExperiment can be created from an existing SOMA", {
  skip_if_not_installed("SummarizedExperiment")
  uri <- file.path(withr::local_tempdir(), "singlecellexperiment")

  # start with scdataset so the soma includes annot matrices
  scdataset <- SOMACollection$new(uri = uri)
  scdataset$from_seurat(pbmc_small)

  soma <- SOMA$new(scdataset$somas$RNA$uri)
  sce <- soma$to_single_cell_experiment(layers = c("counts", "data"))
  expect_s4_class(sce, "SingleCellExperiment")

  expect_setequal(
    names(SingleCellExperiment::reducedDims(sce)),
    SeuratObject::Reductions(pbmc_small)
  )
})
