test_that("SOMASparseNdArray creation", {
  uri <- withr::local_tempdir("sparse-ndarray")
  ndarray <- SOMASparseNdArray$new(uri)
  ndarray$create(arrow::int32(), shape = c(10, 10))

  expect_equal(tiledb::tiledb_object_type(uri), "ARRAY")
  expect_equal(ndarray$dimnames(), c("__dim_0", "__dim_1"))

  expect_equal(ndarray$attrnames(), "data")
  expect_equal(tiledb::datatype(ndarray$attributes()$data), "INT32")

  mat <- create_sparse_matrix_with_int_dims(10, 10)
  ndarray$write(mat)
})
