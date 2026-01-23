# Tests for Carrara Error Handling and Edge Cases ---------------------------

test_that("SOMACollection invalid operations", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()

  # Opening non-existent URI should error
  nonexistent_uri <- file_path(get_carrara_base_uri(), "does_not_exist")
  expect_error(
    SOMACollectionOpen(nonexistent_uri),
    class = "error"
  )

  # Test 2: Adding duplicate member names should error
  SOMACollectionCreate(uri = uri)$close()

  collection <- SOMACollectionOpen(uri, mode = "WRITE")

  # Add first member
  child1_uri <- file.path(uri, "duplicate_name")
  SOMACollectionCreate(child1_uri)$close()

  expect_true("duplicate_name" %in% collection$names())

  # Attempt to add another member with same name should fail
  child2_uri <- file.path(uri, "duplicate_name_v2")
  child2 <- SOMACollectionCreate(child2_uri)

  expect_error(
    collection$add_new_collection(child2, "duplicate_name"),
    class = "error"
  )

  child2$close()
  collection$close()

  # Test 3: Operations on closed objects should error
  collection <- SOMACollectionOpen(uri, mode = "READ")
  collection$close()

  # Attempting to use closed collection should error
  expect_error(
    collection$names(),
    class = "error"
  )
})

test_that("SOMACollection relative URI behavior", {
  skip_if_no_carrara()
  with_carrara_env()

  # Create nested collection structure
  base_uri <- get_carrara_base_uri()
  parent_uri <- file_path(base_uri, generate_unique_id("parent-"))
  child_name <- "child_collection"
  child_uri <- file.path(parent_uri, child_name)

  # Register cleanup for parent (will recursively clean child)
  withr::defer({
    tryCatch(
      {
        grp <- tiledb::tiledb_group(parent_uri)
        tiledb::tiledb_group_close(grp)
        grp <- tiledb::tiledb_group_open(grp, type = "MODIFY_EXCLUSIVE")
        tiledb::tiledb_group_delete(grp = grp, uri = parent_uri, recursive = TRUE)
        tiledb::tiledb_group_close(grp)
      },
      error = function(e) {
        message("Failed to cleanup: ", parent_uri)
      }
    )
  })

  # Create parent and child
  SOMACollectionCreate(parent_uri)$close()
  SOMACollectionCreate(child_uri)$close()

  # Verify child is accessible from parent
  parent <- SOMACollectionOpen(parent_uri, mode = "READ")
  expect_true(child_name %in% parent$names())

  # Verify we can open the child
  child <- parent$get(child_name)
  expect_equivalent(child$soma_type, "SOMACollection")
  child$close()

  parent$close()

  # Verify child can be opened directly via its full URI
  child_direct <- SOMACollectionOpen(child_uri, mode = "READ")
  expect_equivalent(child_direct$soma_type, "SOMACollection")
  child_direct$close()
})

test_that("SOMACollection delete member by name", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri)$close()

  # Add members to the collection
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$add_new_sparse_ndarray(
    "array_to_delete",
    type = arrow::int8(),
    shape = c(10L)
  )
  collection$add_new_collection(
    SOMACollectionCreate(file_path(uri, "collection_to_keep")),
    "collection_to_keep"
  )
  collection$close()

  # Verify both members exist
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_setequal(collection$names(), c("array_to_delete", "collection_to_keep"))
  collection$close()

  # Delete the array member
  collection <- SOMACollectionOpen(uri, mode = "DELETE")
  collection$remove("array_to_delete")
  collection$close()

  # Verify only the collection remains
  collection <- SOMACollectionOpen(uri, mode = "READ")
  expect_equal(collection$names(), "collection_to_keep")
  collection$close()
})

test_that("Invalid nested storage paths are rejected", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()

  # Attempting to create objects with nested storage URIs should fail
  expect_error(
    SOMACollectionCreate(file_path(uri, "s3://bucket/path"))
  )

  expect_error(
    SOMASparseNDArrayCreate(
      uri = file_path(uri, "s3://bucket/array"),
      type = arrow::float32(),
      shape = c(10, 11)
    )
  )
})

test_that("Opening nested child paths directly works", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri)$close()

  # Create nested structure: collection -> c1 -> c1.1 -> dnda1
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$add_new_sparse_ndarray("snda1", type = arrow::int8(), shape = c(10, 11))
  collection$add_new_collection(
    SOMACollectionCreate(file_path(uri, "c1")),
    "c1"
  )
  collection$close()

  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  c1 <- collection$get("c1")
  c1$add_new_collection(
    SOMACollectionCreate(file_path(c1$uri, "c1.1")),
    "c1.1"
  )
  c1$close()
  collection$close()

  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  c1 <- collection$get("c1")
  c1_1 <- c1$get("c1.1")
  c1_1$add_new_dense_ndarray("dnda1", type = arrow::int64(), shape = c(100, 10))
  c1_1$close()
  c1$close()
  collection$close()

  # Verify we can open each nested path directly
  expect_equivalent(SOMACollectionOpen(uri)$soma_type, "SOMACollection")
  expect_equal(
    SOMASparseNDArrayOpen(file_path(uri, "snda1"))$soma_type,
    "SOMASparseNDArray"
  )
  expect_equivalent(
    SOMACollectionOpen(file_path(uri, "c1"))$soma_type,
    "SOMACollection"
  )
  expect_equivalent(
    SOMACollectionOpen(file_path(uri, "c1/c1.1"))$soma_type,
    "SOMACollection"
  )
  expect_equal(
    SOMADenseNDArrayOpen(file_path(uri, "c1/c1.1/dnda1"))$soma_type,
    "SOMADenseNDArray"
  )
})

# Tests for Duplicate Key Handling ------------------------------------------

test_that("SOMACollection add_new_* rejects duplicate key after reopen", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()
  SOMACollectionCreate(uri)$close()

  # Add first member
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$add_new_sparse_ndarray(
    "foo",
    type = arrow::int32(),
    shape = c(5, 5)
  )
  expect_true("foo" %in% collection$names())

  # Attempt to add duplicate in same session
  expect_error(
    collection$add_new_sparse_ndarray(
      "foo",
      type = arrow::int32(),
      shape = c(10, 10)
    ),
    regexp = "Member 'foo' already exists"
  )
  collection$close()

  # Reopen and attempt to add duplicate
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  withr::defer(collection$close())

  expect_true("foo" %in% collection$names())

  expect_error(
    collection$add_new_sparse_ndarray(
      "foo",
      type = arrow::int32(),
      shape = c(10, 10)
    ),
    regexp = "Member 'foo' already exists"
  )
  collection$close()

  # Verify original still there
  collection <- SOMACollectionOpen(uri)
  expect_s3_class(collection$get("foo"), "SOMASparseNDArray")
  expect_equal(
    collection$get("foo")$shape(),
    c(5, 5)
  )
})

test_that("write_soma fails for duplicate keys (same URI)", {
  skip_if_no_carrara()
  with_carrara_env()

  # For Carrara, key must equal URI basename. Attempting to write a second
  # object with the same key means writing to the same URI, which fails because
  # the object already exists.

  uri <- carrara_group_path()
  collection <- SOMACollectionCreate(uri)
  withr::defer(collection$close())

  df1 <- data.frame(a = 1:5, b = letters[1:5])
  df2 <- data.frame(x = 6:10, y = letters[6:10])

  # Write first object (key must match URI basename for Carrara)
  sdf1 <- write_soma(
    df1,
    uri = file_path(uri, "foo"),
    soma_parent = collection,
    key = "foo"
  )
  sdf1$close()
  expect_true("foo" %in% collection$names())

  # Attempt to write another object with same key
  expect_error(
    write_soma(
      df2,
      uri = file_path(uri, "foo"),
      soma_parent = collection,
      key = "foo"
    ),
    regexp = "already exists"
  )

  # Verify original data preserved
  expect_equal(collection$length(), 1L)
  collection$close()

  collection <- SOMACollectionOpen(uri)
  expect_identical(
    collection$get("foo")$read()$concat()$a$as_vector(),
    df1$a
  )
})
