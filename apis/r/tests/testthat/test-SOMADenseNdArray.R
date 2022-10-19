test_that("SOMADenseNdArray creation", {
  # skip_if(TRUE) # temporary
  uri <- withr::local_tempdir("dense-ndarray")

  ndarray <- SOMADenseNdArray$new(uri)
  ndarray$create(arrow::int32(), shape = c(10, 5))

  expect_equal(tiledb::tiledb_object_type(uri), "ARRAY")
  expect_equal(ndarray$dimnames(), c("soma_dim_0", "soma_dim_1"))
  expect_equal(ndarray$attrnames(), "soma_data")
  expect_equal(tiledb::datatype(ndarray$attributes()$soma_data), "INT32")

  mat <- create_dense_matrix_with_int_dims(10, 5)
  ndarray$write(mat)
})
