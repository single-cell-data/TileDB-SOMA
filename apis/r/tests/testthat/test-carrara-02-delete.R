# Tests for Carrara Member Deletion -----------------------------------------

test_that("SOMACollection member deletion", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()

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

  dense <- collection$add_new_dense_ndarray(
    "test_dense",
    type = arrow::float32(),
    shape = c(10L, 10L)
  )
  dense$close()

  child_coll_uri <- file_path(uri, "test_collection")
  child_coll <- SOMACollectionCreate(child_coll_uri)$close()

  # Verify all members exist
  expect_equal(collection$length(), 4)
  expect_setequal(
    collection$names(),
    c("test_dataframe", "test_sparse", "test_dense", "test_collection")
  )

  collection$close()

  # Test deletion of DataFrame
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_no_error(collection$remove("test_dataframe"))
  expect_equal(collection$length(), 3)
  expect_false("test_dataframe" %in% collection$names())
  collection$close()

  # Test deletion of SparseNDArray
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_no_error(collection$remove("test_sparse"))
  expect_equal(collection$length(), 2)
  expect_false("test_sparse" %in% collection$names())
  collection$close()

  # Test deletion of DenseNDArray
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_no_error(collection$remove("test_dense"))
  expect_equal(collection$length(), 1)
  expect_false("test_dense" %in% collection$names())
  collection$close()

  # Test deletion of Collection
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_no_error(collection$remove("test_collection"))
  expect_equal(collection$length(), 0)
  expect_false("test_collection" %in% collection$names())
  collection$close()

  # Test error when trying to remove non-existent member
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_error(
    collection$remove("nonexistent"),
    "does not exist"
  )
  collection$close()
})
