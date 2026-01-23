# Tests for Cloud Error Handling --------------------------------------------

test_that("Error handling for cloud operations", {
  skip_if_no_cloud()

  # Test 1: Opening non-existent URI should error
  nonexistent_uri <- file_path(get_cloud_base_uri(), "does_not_exist_12345")
  expect_error(
    SOMAExperimentOpen(nonexistent_uri),
    class = "error"
  )

  # Test 2: Operations on closed objects should error
  skip_if_not_installed("SeuratObject", .MINIMUM_SEURAT_VERSION("c"))

  pbmc_small <- get_test_seurat_object()
  uri <- cloud_path()
  write_soma(pbmc_small, uri = uri)

  exp <- SOMAExperimentOpen(uri, mode = "READ")
  exp$close()

  # Attempting to use closed experiment should error
  expect_error(
    exp$obs,
    class = "error"
  )
})

# Tests for Duplicate Key Handling --------------------------------------------

test_that("SOMACollection set() rejects duplicate key in same session", {
  skip_if_no_cloud()

  uri <- cloud_path()
  collection <- SOMACollectionCreate(uri)
  withr::defer(collection$close())

  # Create first dataframe
  tbl <- arrow::arrow_table(
    soma_joinid = bit64::as.integer64(0:4),
    value = 1:5
  )

  sdf1_uri <- file_path(uri, "df1")
  sdf1 <- SOMADataFrameCreate(sdf1_uri, tbl$schema, domain = list(soma_joinid = c(0L, 4L)))
  sdf1$write(tbl)
  sdf1$close()

  # Create second dataframe
  tbl2 <- arrow::arrow_table(
    soma_joinid = bit64::as.integer64(0:9),
    value = 6:15
  )
  sdf2_uri <- file_path(uri, "df2")
  sdf2 <- SOMADataFrameCreate(sdf2_uri, tbl2$schema, domain = list(soma_joinid = c(0L, 9L)))
  sdf2$write(tbl2)
  sdf2$close()

  # Set first with key "foo"
  sdf1 <- SOMADataFrameOpen(sdf1_uri)
  collection$set(sdf1, name = "foo")
  expect_true("foo" %in% collection$names())

  # Attempt to set second with same key - should fail
  sdf2 <- SOMADataFrameOpen(sdf2_uri)
  expect_error(
    collection$set(sdf2, name = "foo"),
    regexp = "replacing key 'foo' is unsupported"
  )

  # Attempt add_new_* with same key
  expect_error(
    collection$add_new_sparse_ndarray(key = "foo", type = arrow::int32(), shape = c(10, 10)),
    regexp = "Member 'foo' already exists"
  )

  # Verify original still there
  collection$close()
  collection <- SOMACollectionOpen(uri)
  expect_s3_class(collection$get("foo"), "SOMADataFrame")
  expect_equal(collection$get("foo")$read()$concat()$num_rows, 5)
})

test_that("SOMACollection set() rejects duplicate key after reopen", {
  skip_if_no_cloud()

  uri <- cloud_path()
  collection <- SOMACollectionCreate(uri)

  # Create and add first dataframe
  tbl <- arrow::arrow_table(
    soma_joinid = bit64::as.integer64(0:4),
    value = 1:5
  )
  sdf1_uri <- file_path(uri, "df1")
  sdf1 <- SOMADataFrameCreate(sdf1_uri, tbl$schema, domain = list(soma_joinid = c(0L, 4L)))
  sdf1$write(tbl)
  sdf1$close()

  sdf1 <- SOMADataFrameOpen(sdf1_uri)
  collection$set(sdf1, name = "foo")
  collection$close()

  # Create second dataframe
  tbl2 <- arrow::arrow_table(
    soma_joinid = bit64::as.integer64(0:9),
    value = 6:15
  )
  sdf2_uri <- file_path(uri, "df2")
  sdf2 <- SOMADataFrameCreate(sdf2_uri, tbl2$schema, domain = list(soma_joinid = c(0L, 9L)))
  sdf2$write(tbl2)
  sdf2$close()

  # Reopen and attempt duplicate
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  withr::defer(collection$close())

  expect_true("foo" %in% collection$names())

  sdf2 <- SOMADataFrameOpen(sdf2_uri)
  expect_error(collection$set(sdf2, name = "foo"))

  # Attempt add_new_* with same key after reopen
  expect_error(
    collection$add_new_sparse_ndarray(key = "foo", type = arrow::int32(), shape = c(10, 10)),
    regexp = "Member 'foo' already exists"
  )

  # Verify original still there
  collection$close()
  collection <- SOMACollectionOpen(uri)
  expect_s3_class(collection$get("foo"), "SOMADataFrame")
  expect_equal(collection$get("foo")$read()$concat()$num_rows, 5)
})

test_that("write_soma throws existingKeyWarning for duplicate keys", {
  skip_if_no_cloud()

  uri <- cloud_path()
  collection <- SOMACollectionCreate(uri)
  withr::defer(collection$close())

  df1 <- data.frame(a = 1:5, b = letters[1:5])
  df2 <- data.frame(x = 6:10, y = letters[6:10])

  # Write first object
  sdf1 <- write_soma(df1, uri = file_path(uri, "df1"), soma_parent = collection, key = "foo")
  sdf1$close()
  expect_true("foo" %in% collection$names())

  # Attempt to write another object with same key (verbose = TRUE)
  expect_warning(
    withr::with_options(
      list(verbose = TRUE),
      sdf2 <- write_soma(df2, uri = file_path(uri, "df2"), soma_parent = collection, key = "foo")
    ),
    class = "existingKeyWarning"
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
