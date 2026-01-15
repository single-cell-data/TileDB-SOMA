# Tests for Cloud Object Creation -------------------------------------------

test_that("SOMACollection member operations in cloud", {
  skip_if_no_cloud()
  uri <- cloud_path()

  # Create collection with multiple member types
  SOMACollectionCreate(uri = uri)$close()

  collection <- SOMACollectionOpen(uri, mode = "WRITE")

  # Add different SOMA types as members
  schema <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("value", arrow::int32())
  )
  domain <- list(soma_joinid = c(0L, 100L))
  df <- collection$add_new_dataframe(
    "test_dataframe",
    schema = schema,
    index_column_names = "soma_joinid",
    domain = domain
  )
  df$close()

  sparse <- collection$add_new_sparse_ndarray(
    "test_sparse",
    type = arrow::float32(),
    shape = c(10L, 10L)
  )
  sparse$close()

  # Verify members exist
  expect_equal(collection$length(), 2)
  expect_setequal(
    collection$names(),
    c("test_dataframe", "test_sparse")
  )

  collection$close()

  # Verify members can be opened and have correct types
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_equal(collection$get("test_dataframe")$soma_type, "SOMADataFrame")
  expect_equal(collection$get("test_sparse")$soma_type, "SOMASparseNDArray")
  collection$close()
})
