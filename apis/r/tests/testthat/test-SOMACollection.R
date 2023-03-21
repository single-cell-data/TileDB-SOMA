
test_that("SOMACollection basics", {
  uri <- file.path(withr::local_tempdir(), "new-collection")
  collection <- SOMACollection$new(uri, internal_use_only = "allowed_use")
  expect_equal(collection$uri, uri)

  # Should not exist on disk until created
  expect_false(dir.exists(uri))
  expect_false(collection$exists())

  # Create the collection on disk
  collection$create()

  expect_true(dir.exists(uri))
  expect_match(tiledb::tiledb_object_type(uri), "GROUP")
  expect_true(collection$soma_type == "SOMACollection")
  expect_true(collection$exists())
  expect_equal(collection$length(), 0)

  # Add an element to the collection
  dataframe <- create_and_populate_soma_dataframe(file.path(uri, "sdf"))
  collection$set(dataframe, name = "sdf")

  # Read back the collection
  readback_collection <- SOMACollection$new(uri, internal_use_only = "allowed_use")
  expect_equal(readback_collection$length(), 1)

  readback_dataframe <- readback_collection$get("sdf")
  expect_is(readback_dataframe, "SOMADataFrame")

  # Add the collection to the collection
  collection$add_new_collection(collection, "collection")
  collection2 <- collection$get("collection")
  expect_true(collection2$soma_type == "SOMACollection")
  expect_true(collection2$exists())
  df2 <- collection2$get("sdf")
  ## -- uri differs by "./" so cannot compare R6 objects directly
  expect_equal(readback_dataframe$object[], df2$object[])

  # Add new dataframe to the collection
  collection$add_new_dataframe("new_df", create_arrow_schema(), "foo")
  df3 <- collection$get("new_df")
  expect_true(df3$soma_type == "SOMADataFrame")

  # Add new DenseNdArray to the collection
  collection$add_new_dense_ndarray("nd_d_arr", arrow::int32(), shape = c(10, 5))
  arr <- collection$get("nd_d_arr")
  expect_true(arr$soma_type == "SOMADenseNDArray")

  # Add new SparseNdArray to the collection
  collection$add_new_sparse_ndarray("nd_s_arr", arrow::int32(), shape = c(10, 5))
  arr <- collection$get("nd_s_arr")
  expect_true(arr$soma_type == "SOMASparseNDArray")
})
