# Tests for Carrara Error Handling and Edge Cases --------------------------

test_that("SOMACollection invalid operations", {
  skip_if_no_carrara()
  with_carrara_env()

  uri <- carrara_group_path()

  # Opening non-existent URI should error
  nonexistent_uri <- file_path(get_base_uri(), "does_not_exist")
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
  base_uri <- get_base_uri()
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
