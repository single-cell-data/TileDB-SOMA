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
  collection <- SOMACollectionCreate(uri = uri)$close()

  collection2 <- SOMACollectionOpen(uri)
  # Using equivalence because soma_type returns a named string for collections
  expect_equivalent(collection2$soma_type, "SOMACollection")
  expect_equal(length(collection2$names()), 0)
})

test_that("SOMACollection add_new_* methods", {
  clean_env()
  uri <- create_group_path()

  schema <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("A", arrow::int32())
  )
  domain <- list(soma_joinid = c(0L, 100L))
  type <- arrow::float32()
  shape <- c(99L, 101L)

  children <- paste0("child", 1:7)

  # Create collection
  collection <- SOMACollectionCreate(uri = uri)$close()

  # child1: SOMACollection
  child1_uri <- file.path(uri, "child1")
  child1 <- SOMACollectionCreate(child1_uri)

  # For Carrara URIs children created at nested URIs are automatically added
  # to the parent, so this is a no-op.
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  expect_no_error({
    collection$add_new_collection(child1, "child1")
  })

  # Carrara requires that the member name matches the final segment of the URI
  expect_error(
    collection$add_new_collection(child1, "not_child1"),
    "Member name `not_child1` must match the final segment of the URI"
  )
  collection$close()

  # child2: SOMAExperiment
  child2_uri <- file.path(uri, "child2")
  child2 <- SOMAExperimentCreate(child2_uri)$close()
  child2 <- SOMACollectionOpen(child2_uri)
  expect_true(child2$exists())
  expect_equivalent(child2$soma_type, "SOMAExperiment")

  # child3: SOMAMeasurement
  child3_uri <- file.path(uri, "child3")
  child3 <- SOMAMeasurementCreate(child3_uri)$close()
  child3 <- SOMACollectionOpen(child3_uri)
  expect_true(child3$exists())
  expect_equivalent(child3$soma_type, "SOMAMeasurement")

  # child4: SOMACollection
  child4_uri <- file.path(uri, "child4")
  child4 <- SOMACollectionCreate(child4_uri)$close()
  child4 <- SOMACollectionOpen(child4_uri)
  expect_true(child4$exists())
  expect_equivalent(child4$soma_type, "SOMACollection")

  # child5: SOMADataFrame
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  child5 <- collection$add_new_dataframe(
    "child5d",
    schema = schema,
    index_column_names = "soma_joinid",
    domain = domain
  )
  collection$close()

  # child6: SOMASparseNDArray
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  child6 <- collection$add_new_sparse_ndarray(
    "child6b",
    type = type,
    shape = shape
  )
  collection$close()

  # # child7: SOMADenseNDArray
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  child7 <- collection$add_new_dense_ndarray(
    "child7a",
    type = type,
    shape = shape
  )
  collection$close()

  # # Verify
  # collection <- SOMACollectionOpen(uri, mode = "READ")

  # expect_equal(collection$length(), length(children))
  # expect_setequal(collection$names(), children)

  # # Verify types
  # expect_equal(collection$get("child1")$soma_type, "SOMACollection")
  # expect_equal(collection$get("child2")$soma_type, "SOMAExperiment")
  # expect_equal(collection$get("child3")$soma_type, "SOMAMeasurement")
  # expect_equal(collection$get("child4")$soma_type, "SOMACollection")
  # expect_equal(collection$get("child5")$soma_type, "SOMADataFrame")
  # expect_equal(collection$get("child6")$soma_type, "SOMASparseNDArray")
  # expect_equal(collection$get("child7")$soma_type, "SOMADenseNDArray")

  # # Verify properties for child5, child6, child7
  # c5 <- collection$get("child5")
  # expect_equal(c5$schema(), schema)
  # expect_equal(c5$domain(), domain)

  # c6 <- collection$get("child6")
  # expect_equal(c6$type(), type)
  # expect_equal(c6$shape(), bit64::as.integer64(shape))

  # c7 <- collection$get("child7")
  # expect_equal(c7$type(), type)
  # expect_equal(c7$shape(), bit64::as.integer64(shape))

  # collection$close()
})
