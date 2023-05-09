test_that("Ephemeral Colelction mechanics", {
  uri <- withr::local_tempdir('ephemeral-collection')
  expect_warning(EphemeralCollection$new(uri = uri))
  expect_no_condition(collection <- EphemeralCollection$new())
  expect_true(grepl('^ephemeral-collection:0x[[:digit:]a-f]{6,32}$', collection$uri))
  expect_false(collection$exists())
  expect_s3_class(collection$create(), collection$class())
  expect_true(collection$soma_type == "SOMACollection")
  expect_equal(collection$length(), 0)
  dataframe <- create_and_populate_soma_dataframe(file.path(uri, "sdf"))
  collection$set(dataframe, name = "sdf")
  expect_s3_class(collection$get('sdf'), 'SOMADataFrame')

  # Add the collection to the collection
  expect_error(collection$add_new_collection(collection, "collection"))
  expect_no_condition(collection$set(collection, 'collection'))
  collection2 <- collection$get("collection")
  expect_true(collection2$soma_type == "SOMACollection")
  expect_false(collection2$exists())
  expect_s3_class(df2 <- collection2$get("sdf"), 'SOMADataFrame')

  ## -- uri differs by "./" so cannot compare R6 objects directly
  # expect_equal(readback_dataframe$object[], df2$object[])

  # Add new dataframe to the collection
  expect_error(collection$add_new_dataframe("new_df", create_arrow_schema(), "foo"))
  expect_no_condition(collection$set(
    SOMADataFrameCreate(file.path(uri, 'new_df'), create_arrow_schema(), 'foo'),
    'new_df'
  ))
  expect_s3_class(df3 <- collection$get("new_df"), 'SOMADataFrame')
  expect_true(df3$soma_type == "SOMADataFrame")

  # Add new DenseNdArray to the collection
  expect_error(collection$add_new_dense_ndarray("nd_d_arr", arrow::int32(), shape = c(10, 5)))
  expect_no_condition(collection$set(
    SOMADenseNDArrayCreate(file.path(uri, 'nd_d_arr'), arrow::int32(), c(10, 5)),
    'nd_d_arr'
  ))
  expect_s3_class(arr <- collection$get("nd_d_arr"), 'SOMADenseNDArray')
  expect_true(arr$soma_type == "SOMADenseNDArray")

  # Add new SparseNdArray to the collection
  expect_error(collection$add_new_sparse_ndarray("nd_s_arr", arrow::int32(), shape = c(10, 5)))
  expect_no_condition(collection$set(
    SOMASparseNDArrayCreate(file.path(uri, 'nd_s_arr'), arrow::int32(), c(10, 5)),
    'nd_s_arr'
  ))
  expect_s3_class(arr <- collection$get("nd_s_arr"), 'SOMASparseNDArray')
  expect_true(arr$soma_type == "SOMASparseNDArray")
})
