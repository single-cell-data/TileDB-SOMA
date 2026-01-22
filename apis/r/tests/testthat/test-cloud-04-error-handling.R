# Tests for Cloud Error Handling --------------------------------------------

test_that("Error handling for cloud operations", {
  skip_if_no_cloud()

  # Test 1: Opening non-existent URI should error
  nonexistent_uri <- file_path(get_cloud_base_uri(), "does_not_exist_12345")
  expect_error(
    SOMAExperimentOpen(nonexistent_uri),
    class = "error"
  )

  # Test 2: Operations on closed objects should error
  skip_if_not_installed("SeuratObject", .MINIMUM_SEURAT_VERSION("c"))

  pbmc_small <- get_test_seurat_object()
  uri <- cloud_path()
  write_soma(pbmc_small, uri = uri)

  exp <- SOMAExperimentOpen(uri, mode = "READ")
  exp$close()

  # Attempting to use closed experiment should error
  expect_error(
    exp$obs,
    class = "error"
  )
})
