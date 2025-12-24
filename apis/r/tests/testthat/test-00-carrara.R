# test_that("is_relative_uri checks", {
#   expect_true(is_relative_uri("foo"))
#   expect_true(is_relative_uri("foo/bar"))
#   expect_false(is_relative_uri("/foo/bar"))
#   expect_false(is_relative_uri("file://foo/bar"))
# })

# test_that("sanitize_key checks", {
#   # v3
#   expect_equal(sanitize_key("foo", "tiledbv3"), "foo")
#   # Test that it errors on slash
#   expect_error(sanitize_key("foo/bar", "tiledbv3"), "must not contain slash")

#   # v2
#   expect_equal(sanitize_key("foo", "tiledbv2"), "foo")
#   # Encoded
#   expect_true(sanitize_key("foo/bar", "tiledbv2") != "foo/bar")
# })

PROFILE_NAME <- Sys.getenv("CARRARA_TEST_PROFILE", "carrara")
WORKSPACE_NAME <- Sys.getenv("CARRARA_TEST_WORKSPACE", "TileDB-Inc-Staging")
TEAMSPACE_NAME <- Sys.getenv("CARRARA_TEST_TEAMSPACE", "aaron-dev")
TEST_FOLDER <- Sys.getenv("CARRARA_TEST_FOLDER", "remote_test")
REST_SERVER <- Sys.getenv(
    "TILEDB_REST_SERVER_ADDRESS",
    "https://api.staging.tiledb.io"
)

BASE_URI <- sprintf(
    "tiledb://%s/%s/%s",
    WORKSPACE_NAME,
    TEAMSPACE_NAME,
    TEST_FOLDER
)

clean_env <- function() {
    env_vars <- names(Sys.getenv())
    tiledb_vars <- which(grepl("CARRARA", env_vars) | grepl("TILEDB", env_vars))
    for (var in env_vars[tiledb_vars]) {
        message(sprintf("Unsetting %s", var))
        Sys.unsetenv(var)
    }
}

create_arrow_table <- function(nrows = 10L, factors = FALSE) {
  arrow::arrow_table(
    soma_joinid = seq(bit64::as.integer64(0L), to = nrows - 1L),
    int_column = seq.int(nrows) + 1000L,
    float_column = seq(nrows) + 0.1,
    string_column = as.character(seq.int(nrows) + 1000L),
    schema = create_arrow_schema(FALSE)
  )
}

create_array_path <- function() {
  library(uuid)
  file.path(BASE_URI, paste0(UUIDgenerate(), "-tiledbsoma-r"))
}

create_group_path <- function() {
  library(uuid)
  file.path(BASE_URI, paste0(UUIDgenerate(), "-tiledbsoma-r-group"))
}

test_that("SOMAContext detects data protocol", {
  # Initialize with simple config
  ctx <- SOMAContext$new()

  expect_equal(ctx$get_data_protocol("s3://foo/bar"), "tiledbv2")
  expect_equal(ctx$get_data_protocol("tiledb://foo/bar/baz"), "tiledbv3")

  expect_true(ctx$is_tiledbv2("s3://foo/bar"))
  expect_false(ctx$is_tiledbv3("s3://foo/bar"))
  expect_true(ctx$is_tiledbv3("tiledb://foo/bar/baz"))
  expect_false(ctx$is_tiledbv2("tiledb://foo/bar/baz"))

  # Cloud check
  # ctx_cloud <- SOMAContext$new(c(
  #   "rest.server_address" = "https://api.tiledb.com"
  # ))
  # expect_equal(ctx_cloud$get_data_protocol("tiledb://foo/bar"), "tiledbv2")

  # Carrara check
  # ctx_carrara <- SOMAContext$new(c(
  #   "rest.server_address" = "https://api.something-else.com"
  # ))
  # expect_equal(ctx_carrara$get_data_protocol("tiledb://foo/bar"), "tiledbv3")
})

test_that("SOMADataFrame creation", {
  clean_env()
  uri <- create_array_path()
  table <- create_arrow_table()
  domain <- list(soma_joinid = c(0L, 10L))

  sdf0 <- SOMADataFrameCreate(
    uri = uri,
    schema = table$schema,
    domain = domain
  )
  testthat::expect_no_error({
    sdf0$write(table)
    sdf0$close()
  })
  expect_true(sdf0$exists())

  sdf1 <- SOMADataFrameOpen(uri)
  expect_equal(sdf1$soma_type, "SOMADataFrame")
  expect_equal(sdf1$domain(), domain)
  expect_equal(sdf1$schema(), table$schema)
  expect_equivalent(sdf1$read()$concat(), table)
})

# TODO: Expand to Experiment and Measurements
test_that("SOMACollection creation", {
  clean_env()
  uri <- create_group_path()
  collection <- SOMACollectionCreate(uri = uri)
  collection$close()

  collection2 <- SOMACollectionOpen(uri)
  # Using equivalence because soma_type returns a named string for collections
  expect_equivalent(collection2$soma_type, "SOMACollection")
  expect_equal(length(collection2$names()), 0)
})
