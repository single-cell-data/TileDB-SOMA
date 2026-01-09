# Tests for Carrara URI Enforcement -----------------------------------------

test_that("SOMACollection member name vs URI enforcement", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri = uri)$close()

  # Test 1: Verify error when member name doesn't match URI
  child1_uri <- file_path(uri, "actual_name")
  child1 <- SOMACollectionCreate(child1_uri)

  collection <- SOMACollectionOpen(uri, mode = "WRITE")

  # This should fail because "wrong_name" doesn't match "actual_name"
  expect_error(
    collection$add_new_collection(child1, "wrong_name"),
    "Member name `wrong_name` must match the final segment of the URI"
  )

  collection$close()
  child1$close()

  # Test 2: Verify success when member name matches URI
  child2_uri <- file.path(uri, "correct_name")
  child2 <- SOMACollectionCreate(child2_uri)

  collection <- SOMACollectionOpen(uri, mode = "WRITE")

  # This should succeed because "correct_name" matches URI segment
  expect_no_error(
    collection$add_new_collection(child2, "correct_name")
  )

  expect_true("correct_name" %in% collection$names())
  collection$close()
  child2$close()

  # Test 3: Verify enforcement with special characters in names
  child3_uri <- file.path(uri, "name-with_special.chars")
  child3 <- SOMACollectionCreate(child3_uri)

  collection <- SOMACollectionOpen(uri, mode = "WRITE")

  # Should fail with wrong name
  expect_error(
    collection$add_new_collection(child3, "different_name"),
    "Member name `different_name` must match the final segment of the URI"
  )

  # Should succeed with correct name
  expect_no_error(
    collection$add_new_collection(child3, "name-with_special.chars")
  )

  collection$close()
  child3$close()

  # Test 4: Verify enforcement works for all SOMA types
  collection <- SOMACollectionOpen(uri, mode = "READ")

  # We already tested SOMACollection above, verify it's in the collection
  expect_true("correct_name" %in% collection$names())
  expect_true("name-with_special.chars" %in% collection$names())

  collection$close()
})

test_that("Path separator in member name is rejected", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri = uri)$close()

  # Member names containing path separators should fail with the Carrara because
  # the intermediate path component doesn't exist
  collection <- SOMACollectionOpen(uri, mode = "WRITE")

  schema <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("value", arrow::int32())
  )
  domain <- list(soma_joinid = c(0L, 100L))

  expect_error(
    collection$add_new_dataframe(
      "bad/key",
      schema = schema,
      index_column_names = "soma_joinid",
      domain = domain
    )
  )

  collection$close()
})

test_that("set() is not supported for Carrara collections", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri = uri)$close()

  # Create a child array directly
  child_uri <- file_path(uri, "child_array")
  child <- SOMASparseNDArrayCreate(
    uri = child_uri,
    type = arrow::int32(),
    shape = c(10, 10)
  )

  # set() should fail for Carrara regardless of whether the object exists
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  expect_error(
    collection$set(child, "child_array"),
    class = "unsupportedOperationError"
  )

  # Even with a different name, set() should fail (before name validation)
  expect_error(
    collection$set(child, "different_name"),
    class = "unsupportedOperationError"
  )

  collection$close()
  child$close()
})
