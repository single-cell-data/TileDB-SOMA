
test_that("SOMACollection basics", {
  uri <- file.path(withr::local_tempdir(), "new-collection")

  # Create an empty collection
  collection <- SOMACollectionCreate(uri)
  expect_equal(collection$uri, uri)
  collection$close()

  # Verify the empty collection is accessible and reads back as empty
  collection <- SOMACollectionOpen(uri)
  expect_true(dir.exists(uri))
  expect_match(tiledb::tiledb_object_type(uri), "GROUP")
  expect_true(collection$soma_type == "SOMACollection")
  expect_true(collection$exists())
  expect_equal(collection$length(), 0)
  collection$close()

  # Add a dataframe element to the collection, bypassing add_new_dataframe
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  dataframe <- create_and_populate_soma_dataframe(file.path(uri, "sdf"))
  collection$set(dataframe, name = "sdf")
  collection$close()

  # Read back the collection
  readback_collection <- SOMACollectionOpen(uri)
  expect_equal(readback_collection$length(), 1)
  readback_dataframe <- readback_collection$get("sdf")
  expect_is(readback_dataframe, "SOMADataFrame")
  readback_dataframe$close()
  readback_collection$close()

  # Add a subcollection to the collection
  subcollection <- SOMACollectionCreate(file.path(uri, "subcollection"))$close()

  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$add_new_collection(subcollection, "subcollection")

  subcollection <- collection$get("subcollection")
  subcollection <- SOMACollectionOpen(subcollection$uri)
  expect_true(subcollection$soma_type == "SOMACollection")
  expect_true(subcollection$exists())
  subcollection$close()

  # Add another dataframe to the collection, this time using add_new_dataframe
  collection$add_new_dataframe("new_df", create_arrow_schema(), "foo")$close()
  df3 <- collection$get("new_df")
  df3 <- SOMADataFrameOpen(df3$uri)
  expect_true(df3$soma_type == "SOMADataFrame")
  df3$close()

  # Add new DenseNDArray to the collection
  collection$add_new_dense_ndarray("nd_d_arr", arrow::int32(), shape = c(10, 5))$close()
  arr <- collection$get("nd_d_arr")
  arr <- SOMADenseNDArrayOpen(arr$uri)
  expect_true(arr$soma_type == "SOMADenseNDArray")
  arr$close()

  # Add new SparseNDArray to the collection
  collection$add_new_sparse_ndarray("nd_s_arr", arrow::int32(), shape = c(10, 5))$close()
  arr <- collection$get("nd_s_arr")
  arr <- SOMASparseNDArrayOpen(arr$uri)
  expect_true(arr$soma_type == "SOMASparseNDArray")
  arr$close()

  collection$close()
})

test_that("SOMACollection timestamped ops", {
  # Create a collection @ t0
  uri <- file.path(withr::local_tempdir(), "timestamped-collection")
  collection <- SOMACollectionCreate(uri)
  expect_equal(collection$uri, uri)
  collection$close()
  t0 <- Sys.time()
  Sys.sleep(1.01)

  # add array A with 1 in top-left entry @ t1
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$add_new_sparse_ndarray("A", arrow::int8(), shape = c(2,2))$write(Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2)))
  collection$close()
  t1 <- Sys.time()
  Sys.sleep(1.01)

  # write 1 into bottom-right of A @ t2
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$get("A")$write(Matrix::sparseMatrix(i = 2, j = 2, x = 1, dims = c(2, 2)))
  collection$close()

  # open A via collection with no timestamp => A should reflect the final state
  collection <- SOMACollectionOpen(uri)
  a <- collection$get("A")$read()$sparse_matrix()$concat()
  expect_equal(sum(a), 2)
  collection$close()

  # open A via collection @ t1 => the last write should not be visible
  collection <- SOMACollectionOpen(uri, tiledb_timestamp = t1)
  expect_true("A" %in% collection$names())
  a <- collection$get("A")$read()$sparse_matrix()$concat()
  expect_equal(sum(a), 1)
  collection$close()

  # open collection @ t0 => A should not even be there
  collection <- SOMACollectionOpen(uri, tiledb_timestamp = t0)
  expect_false("A" %in% collection$names())
  expect_error(collection$get("A"))
  collection$close()
})
