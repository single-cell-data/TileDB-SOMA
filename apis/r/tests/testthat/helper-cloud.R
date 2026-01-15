# Cloud Test Configuration and Helpers ----------------------------------

# Skip cloud tests unless explicitly enabled via environment variable
skip_if_no_cloud <- function() {
  testthat::skip_if_not(
    Sys.getenv("SOMA_TEST_CLOUD", "false") == "true",
    "Cloud tests not enabled. Set SOMA_TEST_CLOUD=true to run."
  )
  testthat::skip_if_not(
    nzchar(Sys.getenv("TILEDB_REST_TOKEN")),
    "TILEDB_REST_TOKEN not set; skipping cloud tests."
  )
}

# Read from environment variables or use defaults
get_cloud_config <- function() {
  list(
    namespace = Sys.getenv("CLOUD_TEST_NAMESPACE", "aaronwolen"),
    bucket = Sys.getenv("CLOUD_TEST_BUCKET", "s3://tiledb-aaron/tiledb-cloud/soma-tests")
  )
}

# Build base URI for cloud tests
get_cloud_base_uri <- function() {
  cfg <- get_cloud_config()
  sprintf("tiledb://%s/%s", cfg$namespace, cfg$bucket)
}

# Create a unique cloud uri with automatic cleanup
cloud_path <- function(env = parent.frame()) {
  remote_path(get_cloud_base_uri(), "tiledbsoma-r-test-", cleanup_group, env)
}

# Get a simplified Seurat object for faster testing
get_test_seurat_object <- function() {
  # Load pbmc_small from SeuratObject package
  pbmc_small <- get_data("pbmc_small", package = "SeuratObject")

  # Remove command logs to simplify
  for (cmd in SeuratObject::Command(pbmc_small)) {
    pbmc_small[[cmd]] <- NULL
  }

  # Clear counts and scale.data to reduce size
  pbmc_small[["RNA"]]@counts <- new("matrix")
  pbmc_small[["RNA"]]@scale.data <- new("matrix")

  # Remove existing reductions (we'll add PCA in tests if needed)
  pbmc_small@reductions <- list()

  pbmc_small
}
