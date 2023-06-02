
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

test_that("Platform config and context are respected by add_ methods", {
  uri <- file.path(withr::local_tempdir(), "new-collection")

  # Set params in the config and context
  cfg <- PlatformConfig$new()
  cfg$set("tiledb", "test", "foo", "bar")
  cfg$get("tiledb", "test", "foo")

  ctx <- SOMATileDBContext$new()
  ctx$set("foo", "bar")
  ctx$get("foo")

  # Create an empty collection
  collection <- SOMACollectionCreate(
    uri = uri,
    platform_config = cfg,
    tiledbsoma_ctx = ctx
  )

  # Add a dataframe element to the collection
  tbl <- create_arrow_table()
  sdf1 <- collection$add_new_dataframe("sdf1", tbl$schema, "soma_joinid")
  sdf1$write(tbl)
  collection$close()

  # Verify the config and context params were inherited
  collection$open("READ", internal_use_only = "allowed_use")
  expect_equal(
    collection$get("sdf1")$platform_config$get("tiledb", "test", "foo"),
    "bar"
  )
  expect_equal(
    collection$get("sdf1")$tiledbsoma_ctx$get("foo"),
    "bar"
  )
  collection$close()

  # Method-level config params override instance params
  collection$open("WRITE", internal_use_only = "allowed_use")
  cfg$set("tiledb", "test", "foo", "baz")
  sdf2 <- collection$add_new_dataframe(
    key = "sdf2",
    schema = tbl$schema,
    index_column_names = "soma_joinid",
    platform_config = cfg
  )
  sdf2$write(tbl)

  expect_equal(
    collection$get("sdf2")$platform_config$get("tiledb", "test", "foo"),
    "baz"
  )
  collection$close()
})
