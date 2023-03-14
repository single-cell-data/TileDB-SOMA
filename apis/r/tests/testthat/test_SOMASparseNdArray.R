test_that("SOMASparseNDArray creation", {
  uri <- withr::local_tempdir("sparse-ndarray")
  ndarray <- SOMASparseNDArray$new(uri)
  ndarray$create(arrow::int32(), shape = c(10, 10))

  expect_equal(tiledb::tiledb_object_type(uri), "ARRAY")
  expect_equal(ndarray$dimnames(), c("soma_dim_0", "soma_dim_1"))

  expect_equal(ndarray$attrnames(), "soma_data")
  expect_equal(tiledb::datatype(ndarray$attributes()$soma_data), "INT32")

  mat <- create_sparse_matrix_with_int_dims(10, 10)
  ndarray$write(mat)

  tbl <- ndarray$read_arrow_table(result_order = "COL_MAJOR")
  expect_true(is_arrow_table(tbl))
  expect_equal(tbl$ColumnNames(), c("soma_dim_0", "soma_dim_1", "soma_data"))

  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    ## need to convert to Csparsematrix first to get x values sorted appropriately
    as.numeric(as(mat, "CsparseMatrix")@x)
  )

  mat2 <- ndarray$read_sparse_matrix(repr="T")
  ## not sure why all.equal(mat, mat2) does not pass
  all.equal(as.numeric(mat), as.numeric(mat2))

  # Subset both dims
  tbl <- ndarray$read_arrow_table(
    coords = list(soma_dim_0=0, soma_dim_1=0:2),
    result_order = "COL_MAJOR"
  )
  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    as.numeric(mat[1, 1:3])
  )

  # Subset both dims, unnamed
  tbl <- ndarray$read_arrow_table(
    coords = list(0, 0:2),
    result_order = "COL_MAJOR"
  )
  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    as.numeric(mat[1, 1:3])
  )

  # Validate TileDB array schema
  arr <- tiledb::tiledb_array(uri)
  sch <- tiledb::schema(arr)
  expect_true(tiledb::is.sparse(sch))
  expect_false(tiledb::allows_dups(sch))

  ## shape
  expect_equal(ndarray$shape(), as.integer64(c(10, 10)))

  ## ndim
  expect_equal(ndarray$ndim(), 2L)

  ## nnz
  expect_equal(ndarray$nnz(), 60L)
})

test_that("SOMASparseNDArray creation with duplicates", {
  uri <- withr::local_tempdir("sparse-ndarray")

  set.seed(42)
  D <- data.frame(rows=sample(100, 10, replace=TRUE),
                  cols=sample(100, 10, replace=TRUE),
                  vals=rnorm(10))

  create_write_check <- function(uri, D, allows_dups, do_dup, expected_nnz) {
      ## write from tiledb "for now"
      dom <- tiledb::tiledb_domain(dims = c(tiledb::tiledb_dim("rows", c(1L, 100L), 100L, "INT32"),
                                            tiledb::tiledb_dim("cols", c(1L, 100L), 100L, "INT32")))
      sch <- tiledb::tiledb_array_schema(dom,
                                         attrs=c(tiledb::tiledb_attr("vals", type = "FLOAT64")),
                                         sparse = TRUE,
                                         allows_dups = allows_dups)
      invisible(tiledb::tiledb_array_create(uri, sch))
      arr <- tiledb::tiledb_array(uri)
      if (do_dup)
          arr[] <- rbind(D, D)
      else
          arr[] <- D

      nda <- SOMASparseNDArray$new(uri)
      expect_equal(nda$nnz(), expected_nnz)

      unlink(uri, recursive=TRUE)
  }

  create_write_check(uri, D, FALSE, FALSE, 10)
  create_write_check(uri, D, TRUE, FALSE, 10)
  create_write_check(uri, D, TRUE, TRUE, 20)
})
