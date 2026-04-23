# Shared Remote Test Utilities

s3_tests <- function() {
  envvars <- c("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_S3_BUCKET")
  return(
    requireNamespace("aws.s3", quietly = TRUE) &&
      all(nzchar(Sys.getenv(envvars)))
  )
}

# Generate a unique ID for test assets
generate_unique_id <- function(pattern = "") {
  basename(tempfile(pattern = pattern))
}

# Cleanup a TileDB group with recursive deletion
# Note: tiledb-r requires this specific sequence of operations
cleanup_group <- function(path) {
  tryCatch(
    {
      grp <- tiledb::tiledb_group(path)
      tiledb::tiledb_group_close(grp)
      grp <- tiledb::tiledb_group_open(grp, type = "MODIFY_EXCLUSIVE")
      tiledb::tiledb_group_delete(grp = grp, uri = path, recursive = TRUE)
      tiledb::tiledb_group_close(grp)
    },
    error = function(e) {
      message("Failed to cleanup group: ", path)
    }
  )
}

# Cleanup a TileDB array/directory
cleanup_array <- function(path) {
  tryCatch(
    tiledb::tiledb_vfs_remove_dir(path),
    error = function(e) {
      message("Failed to cleanup array: ", path)
    }
  )
}

# Activate a TileDB profile for the test scope
#
# Sets TILEDB_PROFILE_NAME and unsets TILEDB_REST_TOKEN to ensure the profile's
# credentials are used without interference from any existing token. Additional
# environment variables can be passed via `...`.
#
# @param profile Name of the TileDB profile to activate
# @param ... Additional environment variables to set (name = value pairs)
# @param env Environment for withr::defer scoping
with_tiledb_profile <- function(profile, ..., env = parent.frame()) {
  envvars <- list(
    TILEDB_PROFILE_NAME = profile,
    TILEDB_REST_TOKEN = NA_character_,
    ...
  )
  withr::local_envvar(envvars, .local_envir = env)
}

# Create a unique remote path with automatic cleanup
# @param base_uri Base URI for the remote storage
# @param prefix Pattern prefix for the unique ID
# @param cleanup_fn Function to call for cleanup (cleanup_group or cleanup_array)
# @param env Environment for withr::defer scoping
remote_path <- function(base_uri, prefix, cleanup_fn, env = parent.frame()) {
  path <- file_path(base_uri, generate_unique_id(prefix))
  withr::defer(cleanup_fn(path), envir = env)
  path
}
