# Tests for Seurat I/O with Carrara URIs -------------------------------------

test_that("write_soma Seurat with no error", {
  skip_if_no_carrara()
  skip_if_not_installed("SeuratObject", .MINIMUM_SEURAT_VERSION("c"))
  with_carrara_env()

  # Get test data and simplify for faster testing
  pbmc_small <- get_data("pbmc_small", package = "SeuratObject")
  # Remove command logs to simplify
  for (cmd in SeuratObject::Command(pbmc_small)) {
    pbmc_small[[cmd]] <- NULL
  }

  # Remove scale.data
  pbmc_small[["RNA"]]@counts <- new("matrix")
  pbmc_small[["RNA"]]@scale.data <- new("matrix")

  # Remove reductions
  pbmc_small@reductions <- list()

  uri <- carrara_group_path()

  # Ingest Seurat object to Carrara
  expect_no_error(result_uri <- write_soma(pbmc_small, uri))
})
