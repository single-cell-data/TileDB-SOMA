test_that("SOMACollection basics", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "new-collection")

  # Create an empty collection
  collection <- SOMACollectionCreate(uri)
  expect_equal(collection$uri, uri)
  collection$close()

  # Verify the empty collection is accessible and reads back as empty
  collection <- SOMACollectionOpen(uri)
  expect_true(dir.exists(uri))
  expect_match(
    get_tiledb_object_type(
      collection$uri,
      collection$.__enclos_env__$private$.context$handle
    ),
    "GROUP"
  )
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
  collection$add_new_dataframe(
    "new_df",
    create_arrow_schema(),
    "int_column",
    domain = list(int_column = c(0, 999))
  )$close()
  df3 <- collection$get("new_df")
  df3 <- SOMADataFrameOpen(df3$uri)
  expect_true(df3$soma_type == "SOMADataFrame")
  df3$close()

  # Add new DenseNDArray to the collection
  collection$add_new_dense_ndarray(
    "nd_d_arr",
    arrow::int32(),
    shape = c(10, 5)
  )$close()
  arr <- collection$get("nd_d_arr")
  arr <- SOMADenseNDArrayOpen(arr$uri)
  expect_true(arr$soma_type == "SOMADenseNDArray")
  arr$close()

  # Add new SparseNDArray to the collection
  collection$add_new_sparse_ndarray(
    "nd_s_arr",
    arrow::int32(),
    shape = c(10, 5)
  )$close()
  arr <- collection$get("nd_s_arr")
  arr <- SOMASparseNDArrayOpen(arr$uri)
  expect_true(arr$soma_type == "SOMASparseNDArray")
  arr$close()

  collection$close()
})

test_that("SOMACollection timestamped ops", {
  skip_if(!extended_tests())
  # Create a collection @ t0
  uri <- tempfile(pattern = "timestamped-collection")
  collection <- SOMACollectionCreate(uri)
  expect_equal(collection$uri, uri)
  collection$close()
  t0 <- Sys.time()
  Sys.sleep(1.01)

  # add array A with 1 in top-left entry @ t1
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$add_new_sparse_ndarray(
    "A",
    arrow::int8(),
    shape = c(2, 2)
  )$write(Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2)))
  collection$close()
  t1 <- Sys.time()
  Sys.sleep(1.01)

  # write 1 into bottom-right of A @ t2
  collection <- SOMACollectionOpen(uri, mode = "WRITE")
  collection$get("A")$write(Matrix::sparseMatrix(
    i = 2,
    j = 2,
    x = 1,
    dims = c(2, 2)
  ))
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
  ## FIXME(once groups done via C++)  expect_equal(sum(a), 1)
  collection$close()

  # open collection @ t0 => A should not even be there
  collection <- SOMACollectionOpen(uri, tiledb_timestamp = t0)
  expect_false("A" %in% collection$names())
  expect_error(collection$get("A"))
  collection$close()
})

test_that("Platform config and context are respected by add_ methods", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "new-collection")

  # Set params in the config and context
  cfg <- PlatformConfig$new()
  cfg$set("tiledb", "test", "int_column", "float_column")
  cfg$get("tiledb", "test", "int_column")

  ctx <- SOMAContext$new(config = c(test_key = "test_value"))

  # Create an empty collection
  collection <- SOMACollectionCreate(
    uri = uri,
    platform_config = cfg,
    context = ctx
  )
  on.exit(collection$close(), add = TRUE, after = FALSE)

  # Add a dataframe element to the collection
  tbl <- create_arrow_table()
  sdf1 <- collection$add_new_dataframe(
    "sdf1",
    tbl$schema,
    "soma_joinid",
    domain = list(soma_joinid = c(0, 999))
  )
  sdf1$write(tbl)

  # Verify the config and context params were inherited
  collection <- collection$reopen("READ")
  expect_equal(
    collection$get("sdf1")$platform_config$get("tiledb", "test", "int_column"),
    "float_column"
  )
  output_config <- collection$get("sdf1")$context$get_config()
  expect_type(output_value <- output_config["test_key"], "character")
  expect_equal(unname(output_value), "test_value")

  # Method-level config params override instance params
  collection <- collection$reopen("WRITE")
  cfg$set("tiledb", "test", "int_column", "string_column")
  sdf2 <- collection$add_new_dataframe(
    key = "sdf2",
    schema = tbl$schema,
    index_column_names = "soma_joinid",
    domain = list(soma_joinid = c(0, 999)),
    platform_config = cfg
  )
  sdf2$write(tbl)

  expect_equal(
    collection$get("sdf2")$platform_config$get("tiledb", "test", "int_column"),
    "string_column"
  )
})
