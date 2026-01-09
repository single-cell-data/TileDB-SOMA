# Helper functions for Carrara tests --------------------------------------

# Replacement helper that matches create_arrow_schema() and respects nullability
create_arrow_table <- function(nrows = 10L, factors = FALSE) {
  arrow::arrow_table(
    soma_joinid = seq(bit64::as.integer64(0L), to = nrows - 1L),
    int_column = seq.int(nrows) + 1000L,
    float_column = seq(nrows) + 0.1,
    string_column = as.character(seq.int(nrows) + 1000L),
    schema = create_arrow_schema(FALSE)
  )
}

# Tests for Carrara Object Creation -----------------------------------------

test_that("SOMAContext detects data protocol", {
  skip_if_no_carrara()
  ctx <- SOMAContext$new()

  expect_equal(ctx$get_data_protocol("s3://foo/bar"), "tiledbv2")
  expect_equal(ctx$get_data_protocol("tiledb://foo/bar/baz"), "tiledbv3")

  expect_true(ctx$is_tiledbv2("s3://foo/bar"))
  expect_false(ctx$is_tiledbv3("s3://foo/bar"))
  expect_true(ctx$is_tiledbv3("tiledb://foo/bar/baz"))
  expect_false(ctx$is_tiledbv2("tiledb://foo/bar/baz"))
})

test_that("SOMADataFrame creation", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_array_path()
  tbl <- create_arrow_table()
  domain <- list(soma_joinid = c(0L, 10L))

  sdf0 <- SOMADataFrameCreate(
    uri = uri,
    schema = tbl$schema,
    domain = domain
  )

  expect_no_error({
    sdf0$write(tbl)
    sdf0$close()
  })
  expect_true(sdf0$exists())

  sdf1 <- SOMADataFrameOpen(uri)
  expect_equal(sdf1$soma_type, "SOMADataFrame")
  expect_equal(sdf1$domain(), domain)
  expect_equal(sdf1$schema(), tbl$schema)

  expect_equivalent(sdf1$read()$concat(), tbl)
  sdf1$close()
})

test_that("SOMACollection creation", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri = uri)$close()

  collection2 <- SOMACollectionOpen(uri)
  # Testing equivalence because soma_type returns a named string for collections
  expect_equivalent(collection2$soma_type, "SOMACollection")
  expect_equal(length(collection2$names()), 0)
  collection2$close()
})

test_that("SOMACollection add_new_* methods", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()

  schema <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("A", arrow::int32())
  )
  domain <- list(soma_joinid = c(0L, 100L))
  type <- arrow::float32()
  shape <- c(99L, 101L)

  children <- c(
    "child1",
    "child2",
    "child3",
    "child4",
    "child5",
    "child6",
    "child7"
  )

  # Create collection
  SOMACollectionCreate(uri = uri)$close()

  # child1: SOMACollection
  child1_uri <- file_path(uri, "child1")
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
  child1$close()

  # child2: SOMAExperiment
  child2_uri <- file.path(uri, "child2")
  SOMAExperimentCreate(child2_uri)$close()
  child2 <- SOMACollectionOpen(child2_uri)
  expect_true(child2$exists())
  expect_equivalent(child2$soma_type, "SOMAExperiment")
  child2$close()

  # child3: SOMAMeasurement
  child3_uri <- file.path(uri, "child3")
  SOMAMeasurementCreate(child3_uri)$close()
  child3 <- SOMACollectionOpen(child3_uri)
  expect_true(child3$exists())
  expect_equivalent(child3$soma_type, "SOMAMeasurement")
  child3$close()

  # child4: SOMACollection
  child4_uri <- file.path(uri, "child4")
  SOMACollectionCreate(child4_uri)$close()
  child4 <- SOMACollectionOpen(child4_uri)
  expect_true(child4$exists())
  expect_equivalent(child4$soma_type, "SOMACollection")
  child4$close()

  # child5: SOMADataFrame
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  child5 <- collection$add_new_dataframe(
    "child5",
    schema = schema,
    index_column_names = "soma_joinid",
    domain = domain
  )
  child5$close()
  collection$close()

  # child6: SOMASparseNDArray
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  child6 <- collection$add_new_sparse_ndarray(
    "child6",
    type = type,
    shape = shape
  )
  child6$close()
  collection$close()

  # child7: SOMADenseNDArray
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  child7 <- collection$add_new_dense_ndarray(
    "child7",
    type = type,
    shape = shape
  )
  child7$close()
  collection$close()

  # Verify all children were added
  collection <- SOMACollectionOpen(uri, mode = "READ")

  expect_equal(collection$length(), length(children))
  expect_setequal(collection$names(), children)

  # Verify types on retrieval
  expect_equivalent(collection$get("child1")$soma_type, "SOMACollection")
  expect_equivalent(collection$get("child2")$soma_type, "SOMAExperiment")
  expect_equivalent(collection$get("child3")$soma_type, "SOMAMeasurement")
  expect_equivalent(collection$get("child4")$soma_type, "SOMACollection")
  expect_equivalent(collection$get("child5")$soma_type, "SOMADataFrame")
  expect_equivalent(collection$get("child6")$soma_type, "SOMASparseNDArray")
  expect_equivalent(collection$get("child7")$soma_type, "SOMADenseNDArray")

  # Verify properties for child5, child6, child7
  c5 <- collection$get("child5")
  expect_equal(c5$schema(), schema)
  expect_equal(c5$domain(), domain)
  c5$close()

  c6 <- collection$get("child6")
  expect_equal(c6$shape(), bit64::as.integer64(shape))
  c6$close()

  c7 <- collection$get("child7")
  expect_equal(c7$shape(), bit64::as.integer64(shape))
  c7$close()

  collection$close()

  # Verify that all child URIs are relative
  grp <- tiledb::tiledb_group(uri)
  for (i in seq_len(tiledb::tiledb_group_member_count(grp)) - 1L) {
    child_name <- tiledb::tiledb_group_member(grp, i)[3]
    expect_true(tiledb::tiledb_group_is_relative(grp, child_name))
  }
  tiledb::tiledb_group_close(grp)
})

test_that("SOMACollection path encoding validation", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri = uri)$close()

  # Test names with special characters
  test_names <- c(
    "my collection",           # spaces
    "data_αβγ",                # unicode
    "test-data_v1",            # hyphens and underscores
    "name.with.dots",          # dots
    "CamelCaseName"            # mixed case
  )

  collection <- SOMACollectionOpen(uri, mode = "WRITE")

  # Create child collections with special character names
  for (name in test_names) {
    child_uri <- file_path(uri, name)
    SOMACollectionCreate(child_uri)$close()
  }

  collection$close()

  # Verify all collections were created and names preserved
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_equal(collection$length(), length(test_names))
  expect_setequal(collection$names(), test_names)

  # Verify each collection can be opened and has correct type
  for (name in test_names) {
    child <- collection$get(name)
    expect_equivalent(child$soma_type, "SOMACollection")
    child$close()
  }

  collection$close()

  # Test round-trip: create DataFrame with unicode name
  unicode_df_uri <- file.path(uri, "dataframe_αβγ")
  schema <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("value", arrow::int32())
  )
  domain <- list(soma_joinid = c(0L, 100L))

  SOMADataFrameCreate(
    uri = unicode_df_uri,
    schema = schema,
    domain = domain
  )$close()

  # Verify DataFrame with unicode name
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_true("dataframe_αβγ" %in% collection$names())

  df <- collection$get("dataframe_αβγ")
  expect_equivalent(df$soma_type, "SOMADataFrame")
  df$close()

  collection$close()
})

test_that("SOMACollection reflects carrara's implicit member addition", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri)$close()

  # Create an NDArray child directly (not using add_new_*)
  child_uri <- file_path(uri, "sparse_ndarray_A")
  SOMASparseNDArrayCreate(
    uri = child_uri,
    type = arrow::int32(),
    shape = c(10, 10)
  )$close()

  # Open collection and verify implicit addition
  collection <- SOMACollectionOpen(uri)
  expect_true("sparse_ndarray_A" %in% collection$names())
  expect_true(collection$get("sparse_ndarray_A")$exists())
  collection$close()

  # Carrara does not support the set() operation
  array_A <- SOMASparseNDArrayOpen(child_uri)
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  expect_error(
    collection$set(array_A, "sparse_ndarray_A"),
    class = "unsupportedOperationError"
  )
  collection$close()
  array_A$close()
})
