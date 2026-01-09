# Carrara Test Configuration and Helpers -----------------------------------

# Generate a unique ID for test assets
generate_unique_id <- function(pattern = "") {
  basename(tempfile(pattern = pattern))
}

# Skip carrara tests unless explicitly enabled via environment variable
skip_if_no_carrara <- function() {
  if (Sys.getenv("SOMA_TEST_CARRARA", "false") != "true") {
    skip("Carrara tests not enabled. Set SOMA_TEST_CARRARA=true to run.")
  }
}

# Read from environment variables or use defaults
get_carrara_config <- function() {
  list(
    profile = Sys.getenv("CARRARA_TEST_PROFILE", "carrara"),
    workspace = Sys.getenv("CARRARA_TEST_WORKSPACE", "TileDB-Inc-Staging"),
    teamspace = Sys.getenv("CARRARA_TEST_TEAMSPACE", "aaron-dev"),
    folder = Sys.getenv("CARRARA_TEST_FOLDER", "remote_test"),
    rest_server = Sys.getenv(
      "TILEDB_REST_SERVER_ADDRESS",
      "https://api.staging.tiledb.io"
    )
  )
}

# Set carrara-related environment variables and unset existing TILEDB_REST_TOKEN
with_carrara_env <- function(env = parent.frame()) {
    withr::local_envvar(list(
      CARRARA_PROFILE = get_carrara_config()$profile,
      TILEDB_REST_SERVER_ADDRESS = get_carrara_config()$rest_server,
      TILEDB_REST_TOKEN = NA_character_
  ), .local_envir = env)
}

# Build base URI for carrara tests
get_base_uri <- function() {
  cfg <- get_carrara_config()
  sprintf("tiledb://%s/%s/%s", cfg$workspace, cfg$teamspace, cfg$folder)
}

# Create a unique carrara array path with automatic cleanup
carrara_array_path <- function(env = parent.frame()) {

  path <- file_path(
    get_base_uri(),
    generate_unique_id("tiledbsoma-r-")
  )

  # Register cleanup - delete array after test completes
  withr::defer(
    {
      tryCatch(
        tiledb::tiledb_vfs_remove_dir(path),
        error = function(e) {
          message("Failed to cleanup carrara array: ", path)
        }
      )
    },
    envir = env
  )

  path
}

# Create a unique carrara group path with automatic cleanup
carrara_group_path <- function(env = parent.frame()) {
  path <- file_path(
    get_base_uri(),
    generate_unique_id("tiledbsoma-r-group-")
  )

  # Recursively delete group after test completes
  # Note: tiledb-r requires this specific sequence of operations
  withr::defer(
    {
      tryCatch(
        {
          grp <- tiledb::tiledb_group(path)
          tiledb::tiledb_group_close(grp)
          grp <- tiledb::tiledb_group_open(grp, type = "MODIFY_EXCLUSIVE")
          tiledb::tiledb_group_delete(grp = grp, uri = path, recursive = TRUE)
          tiledb::tiledb_group_close(grp)
        },
        error = function(e) {
          message("Failed to cleanup carrara group: ", path)
        }
      )
    },
    envir = env
  )

  path
}
