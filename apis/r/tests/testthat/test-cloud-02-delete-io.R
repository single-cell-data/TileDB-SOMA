# Tests for Cloud Member Deletion -------------------------------------------

test_that("DELETE mode member removal in cloud", {
  skip_if_no_cloud()
  uri <- cloud_path()

  # Create collection with members
  SOMACollectionCreate(uri = uri)$close()

  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$add_new_sparse_ndarray(
    "array_to_delete",
    type = arrow::int8(),
    shape = c(10L)
  )$close()
  collection$add_new_sparse_ndarray(
    "array_to_keep",
    type = arrow::int8(),
    shape = c(10L)
  )$close()
  collection$close()

  # Verify both members exist
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_setequal(collection$names(), c("array_to_delete", "array_to_keep"))
  collection$close()

  # Delete one member
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_no_error(collection$remove("array_to_delete"))
  collection$close()

  # Verify only one member remains
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_equal(collection$names(), "array_to_keep")
  expect_equal(collection$length(), 1)
  collection$close()

  # Test error when trying to remove non-existent member
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  expect_error(
    collection$remove("nonexistent"),
    "does not exist"
  )
  collection$close()
})
