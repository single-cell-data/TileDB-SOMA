test_that("SOMADenseNDArray creation", {
  uri <- withr::local_tempdir("dense-ndarray")

  ndarray <- SOMADenseNDArray$new(uri, internal_use_only = "allowed_use")
  ndarray$create(arrow::int32(), shape = c(10, 5))

  expect_equal(tiledb::tiledb_object_type(uri), "ARRAY")
  expect_equal(ndarray$dimnames(), c("soma_dim_0", "soma_dim_1"))
  expect_equal(ndarray$attrnames(), "soma_data")
  expect_equal(tiledb::datatype(ndarray$attributes()$soma_data), "INT32")

  mat <- create_dense_matrix_with_int_dims(10, 5)
  ndarray$write(mat)

  # Read result in column-major order to match R matrix layout
  tbl <- ndarray$read_arrow_table(result_order = "COL_MAJOR")
  expect_true(is_arrow_table(tbl))
  expect_equal(tbl$ColumnNames(), c("soma_dim_0", "soma_dim_1", "soma_data"))

  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    as.numeric(mat)
  )

  expect_equal(ndarray$read_dense_matrix(), mat)

  # Subset the array on both dimensions
  tbl <- ndarray$read_arrow_table(
    coords = list(soma_dim_0=0:3, soma_dim_1=0:2),
    result_order = "COL_MAJOR"
  )
  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    as.numeric(mat[1:4, 1:3])
  )

  # Subset the array on both dimensions, unnamed list
  tbl <- ndarray$read_arrow_table(
    coords = list(0:3, 0:2),
    result_order = "COL_MAJOR"
  )
  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    as.numeric(mat[1:4, 1:3])
  )


  # Subset the array on the second dimension
  tbl <- ndarray$read_arrow_table(
    coords = list(soma_dim_1 = bit64::as.integer64(0:2)),
    result_order = "COL_MAJOR"
  )
  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    as.numeric(mat[, 1:3])
  )

  # Validating coords format
  expect_error(
    ndarray$read_arrow_table(coords = list(cbind(0, 1))),
    "must be a list of vectors"
  )

  # Validate TileDB array schema
  arr <- tiledb::tiledb_array(uri)
  sch <- tiledb::schema(arr)
  expect_false(tiledb::is.sparse(sch))

  ## shape
  expect_equal(ndarray$shape(), as.integer64(c(10, 5)))

  ## ndim
  expect_equal(ndarray$ndim(), 2L)

})
