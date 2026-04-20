# Cloud Test Configuration and Helpers ----------------------------------

# Skip cloud tests unless explicitly enabled via environment variable
skip_if_no_cloud <- function() {
  testthat::skip_if_not(
    Sys.getenv("SOMA_TEST_CLOUD", "false") == "true",
    "Cloud tests not enabled. Set SOMA_TEST_CLOUD=true to run."
  )
  testthat::skip_if_not(
    nzchar(Sys.getenv("CLOUD_TEST_BUCKET")),
    "CLOUD_TEST_BUCKET not set; skipping cloud tests."
  )
}

# Read from environment variables or use defaults
get_cloud_config <- function() {
  list(
    profile = Sys.getenv("CLOUD_TEST_PROFILE", "default"),
    namespace = Sys.getenv("CLOUD_TEST_NAMESPACE", "TileDB-Inc"),
    bucket = Sys.getenv("CLOUD_TEST_BUCKET")
  )
}

# Activate the cloud TileDB profile for the test scope
with_cloud_env <- function(env = parent.frame()) {
  with_tiledb_profile(get_cloud_config()$profile, env = env)
}

# Build base URI for cloud tests
get_cloud_base_uri <- function() {
  cfg <- get_cloud_config()
  sprintf("tiledb://%s/%s", cfg$namespace, cfg$bucket)
}

# Create a unique cloud array path with automatic cleanup.
#
# NOTE: tiledb-r does not provide a method for deleting and
# deregistering arrays from TileDB Cloud. The cleanup function
# only removes the array from S3 via VFS. Registered assets
# must be cleaned up manually in the TileDB Cloud UI.
cloud_array_path <- function(env = parent.frame()) {
  remote_path(get_cloud_base_uri(), "tiledbsoma-r-", cleanup_array, env)
}

# Create a unique cloud group path with automatic cleanup
cloud_group_path <- function(env = parent.frame()) {
  remote_path(get_cloud_base_uri(), "tiledbsoma-r-group-", cleanup_group, env)
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
