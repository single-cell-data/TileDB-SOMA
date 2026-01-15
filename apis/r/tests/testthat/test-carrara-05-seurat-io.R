# Tests for Seurat I/O with Carrara URIs -------------------------------------

test_that("write_soma Seurat with no error", {
  skip_if_no_carrara()
  skip_if_not_installed("SeuratObject", .MINIMUM_SEURAT_VERSION("c"))
  with_carrara_env()

  pbmc_small <- get_test_seurat_object()

  uri <- carrara_group_path()

  # Ingest Seurat object to Carrara
  expect_no_error(result_uri <- write_soma(pbmc_small, uri))
})
