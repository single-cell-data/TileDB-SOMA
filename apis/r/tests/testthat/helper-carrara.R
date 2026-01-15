# Carrara Test Configuration and Helpers --------------------------------

# Skip carrara tests unless explicitly enabled via environment variable
skip_if_no_carrara <- function() {
  testthat::skip_if_not(
    Sys.getenv("SOMA_TEST_CARRARA", "false") == "true",
    "Carrara tests not enabled. Set SOMA_TEST_CARRARA=true to run."
  )
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
    TILEDB_PROFILE_NAME = get_carrara_config()$profile,
    TILEDB_REST_SERVER_ADDRESS = get_carrara_config()$rest_server,
    TILEDB_REST_TOKEN = NA_character_
  ), .local_envir = env)
}

# Build base URI for carrara tests
get_carrara_base_uri <- function() {
  cfg <- get_carrara_config()
  sprintf("tiledb://%s/%s/%s", cfg$workspace, cfg$teamspace, cfg$folder)
}

# Create a unique carrara array path with automatic cleanup
carrara_array_path <- function(env = parent.frame()) {
  remote_path(get_carrara_base_uri(), "tiledbsoma-r-", cleanup_array, env)
}

# Create a unique carrara group path with automatic cleanup
carrara_group_path <- function(env = parent.frame()) {
  remote_path(get_carrara_base_uri(), "tiledbsoma-r-group-", cleanup_group, env)
}
