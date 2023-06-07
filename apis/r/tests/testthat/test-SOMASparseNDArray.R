test_that("SOMASparseNDArray creation", {
  uri <- withr::local_tempdir("sparse-ndarray")
  ndarray <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(10, 10))

  expect_equal(tiledb::tiledb_object_type(uri), "ARRAY")
  expect_equal(ndarray$dimnames(), c("soma_dim_0", "soma_dim_1"))

  expect_equal(ndarray$attrnames(), "soma_data")
  expect_equal(tiledb::datatype(ndarray$attributes()$soma_data), "INT32")

  mat <- create_sparse_matrix_with_int_dims(10, 10)
  vals <- as.vector(t(as.matrix(mat)))
  vals <- vals[ vals != 0 ] # needed below for comparison
  ndarray$write(mat)
  ndarray$close()

  ndarray <- SOMASparseNDArrayOpen(uri)

  tbl <- ndarray$read(result_order = "COL_MAJOR")$tables()$concat()
  expect_true(is_arrow_table(tbl))
  expect_equal(tbl$ColumnNames(), c("soma_dim_0", "soma_dim_1", "soma_data"))
  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    ## need to convert to Csparsematrix first to get x values sorted appropriately
    as.numeric(as(mat, "CsparseMatrix")@x)
  )

  ## Subset both dims
  tbl <- ndarray$read(
    coords = list(soma_dim_0=0, soma_dim_1=0:2),
    result_order = "COL_MAJOR"
  )$tables()$concat()
  expect_identical(
    as.numeric(tbl$GetColumnByName("soma_data")),
    as.numeric(mat[1, 1:3])
  )

  ## Subset both dims, unnamed
  tbl <- ndarray$read(
    coords = list(0, 0:2),
    result_order = "COL_MAJOR"
  )$tables()$concat()
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

  ## nnz as free function
  expect_equal(nnz(uri), 60L)
  ## nzz with config, expected breakge as 'bad key' used
  expect_error(nnz(uri, c(sm.encryption_key="Nope", sm.encryption_type="AES_256_GCM")))
  ## shape as free function
  expect_equal(shape(uri), c(10,10))
  ## shape with config, expected breakge as 'bad key' used
  expect_error(shape(uri, c(sm.encryption_key="Nope", sm.encryption_type="AES_256_GCM")))

  ndarray$close()

})

test_that("SOMASparseNDArray read", {
  uri <- withr::local_tempdir("sparse-ndarray")
  ndarray <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(10, 10))

  # For this test, write 9x9 data into 10x10 array. Leaving the last row & column
  # empty touches corner cases with setting dims() correctly
  mat <- create_sparse_matrix_with_int_dims(9, 9, repr = "T")
  ndarray$write(mat)
  expect_equal(as.numeric(ndarray$shape()), c(10, 10))
  ndarray$close()
  
  ndarray <- SOMASparseNDArrayOpen(uri)

  # read sparse matrix directly
  mat2 <- ndarray$read()$sparse_matrix(zero_based = F)$concat()
  expect_s4_class(mat2, "sparseMatrix")
  expect_s4_class(mat2, "dgTMatrix")
  expect_equal(dim(mat2), c(10, 10))
  expect_equal(nrow(mat2), 10)
  expect_equal(ncol(mat2), 10)
  all.equal(mat, mat2) 
  all.equal(mat[1:9, 1:9], mat2[1:9, 1:9])
  
  # test coerced matrix: dgTMatrix
  mat2 <- ndarray$read()$sparse_matrix(zero_based = F)$concat()
  mat_coerced <- as(ndarray$read(), "TsparseMatrix")
  expect_s4_class(mat_coerced, "dgTMatrix")
  all.equal(mat2, mat_coerced) 
  
  # test coerced matrix: dgCMatrix
  mat2 <- ndarray$read()$sparse_matrix(zero_based = F)$concat()
  mat_coerced <- as(ndarray$read(), "CsparseMatrix")
  expect_s4_class(mat_coerced, "dgCMatrix")
  all.equal(as(mat2, "CsparseMatrix"), mat_coerced) 
  
  # test coerced matrix: dgTMatrix
  mat2 <- ndarray$read()$sparse_matrix(zero_based = F)$concat()
  mat_coerced <- as(ndarray$read(), "RsparseMatrix")
  expect_s4_class(mat_coerced, "dgRMatrix")
  all.equal(as(mat2, "RsparseMatrix"), mat_coerced) 
  

  # test zero-based Matrix view
  mat2 <- ndarray$read()$sparse_matrix(zero_based=T)$concat()
  expect_true(inherits(mat2, "matrixZeroBasedView"))
  expect_s4_class(mat2$get_one_based_matrix(), "dgTMatrix")
  expect_equal(mat2$dim(), c(10, 10))
  expect_equal(mat2$nrow(), 10)
  expect_equal(mat2$ncol(), 10)
  all.equal(mat, mat2$get_one_based_matrix)
  all.equal(mat, mat2$take(0:8,0:8)$get_one_based_matrix())
  
  # test zero-based view and convert to sparse matrix
  mat2 <- ndarray$read()$sparse_matrix(zero_based = T)$concat()
  expect_true(inherits(mat2, "matrixZeroBasedView"))
  expect_s4_class(mat2$get_one_based_matrix(), "sparseMatrix")
  expect_equal(mat2$dim(), c(10, 10))
  expect_equal(mat2$nrow(), 10)
  expect_equal(mat2$ncol(), 10)
  all.equal(mat, mat2) 
  all.equal(mat[1:9, 1:9], mat2$take(0:8, 0:8)$get_one_based_matrix())

  # test with iterated reader
  iterator <- ndarray$read()$sparse_matrix(zero_based = T)
  mat2 <- iterator$read_next()
  expect_true(inherits(mat2, "matrixZeroBasedView"))
  expect_s4_class(mat2$get_one_based_matrix(), "sparseMatrix")
  expect_equal(mat2$dim(), c(10, 10))
  expect_equal(mat2$nrow(), 10)
  expect_equal(mat2$ncol(), 10)
  all.equal(mat, mat2$get_one_based_matrix)
  all.equal(mat, mat2$take(0:8,0:8)$get_one_based_matrix())
  
  # test arrow table reader
  df <- data.frame(soma_dim_0 = mat@i, soma_dim_1 = mat@j, soma_data = mat@x)
  tbl <- arrow::arrow_table(df[order(df$soma_dim_0, df$soma_dim_1),])
  tbl2 <- ndarray$read(result_order = "ROW_MAJOR")$tables()$concat()
  all.equal(tbl, tbl2)
  
  # test arrow table coercion
  df <- data.frame(soma_dim_0 = mat@i, soma_dim_1 = mat@j, soma_data = mat@x)
  tbl <- arrow::arrow_table(df[order(df$soma_dim_0, df$soma_dim_1),])
  tbl2 <- arrow::as_arrow_table(ndarray$read(result_order = "ROW_MAJOR"))
  all.equal(tbl, tbl2)
  
  # test data.frame coercion
  df <- data.frame(soma_dim_0 = mat@i, soma_dim_1 = mat@j, soma_data = mat@x)
  df <- df[order(df$soma_dim_0, df$soma_dim_1),]
  df2 <- as.data.frame(ndarray$read(result_order = "ROW_MAJOR"))
  expect_true(inherits(df2, "data.frame"))
  all.equal(df, df2)
  
  ndarray$close()
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

      nda <- SOMASparseNDArray$new(uri, internal_use_only = "allowed_use")
      expect_equal(nda$nnz(), expected_nnz)

      unlink(uri, recursive=TRUE)
  }

  create_write_check(uri, D, FALSE, FALSE, 10)
  create_write_check(uri, D, TRUE, FALSE, 10)
  create_write_check(uri, D, TRUE, TRUE, 20)
})

test_that("platform_config is respected", {
  uri <- withr::local_tempdir("soma-sparse-nd-array")

  # Set tiledb create options
  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'sparse_nd_array_dim_zstd_level', 9)
  cfg$set('tiledb', 'create', 'capacity', 8000)
  cfg$set('tiledb', 'create', 'tile_order', 'COL_MAJOR')
  cfg$set('tiledb', 'create', 'cell_order', 'UNORDERED')
  cfg$set('tiledb', 'create', 'offsets_filters', list("RLE"))
  cfg$set('tiledb', 'create', 'validity_filters', list("RLE", "NONE"))
  cfg$set('tiledb', 'create', 'dims', list(
    soma_dim_0 = list(
      filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=8), "NONE")
      # TODO: test setting/checking tile extent, once shapes/domain-maxes are made programmable.
      # At present we get:
      #
      #   Error: Tile extent check failed; domain max expanded to multiple of tile extent exceeds
      #   max value representable by domain type
      #
      # tile = 999
    ),
    soma_dim_1 = list(
      filters = list("RLE")
      # TODO: test setting/checking tile extent, once shapes/domain-maxes are made programmable.
      # At present we get:
      #
      #   Error: Tile extent check failed; domain max expanded to multiple of tile extent exceeds
      #   max value representable by domain type
      #
      # tile = 999
    )
  ))
  cfg$set('tiledb', 'create', 'attrs', list(
    soma_data = list(
      filters = list("BITSHUFFLE", list(name="ZSTD", COMPRESSION_LEVEL=9))
    )
  ))

  # Create the SOMASparseNDArray
  snda <- SOMASparseNDArrayCreate(uri=uri, type=arrow::int32(), shape=c(100,100), platform_config = cfg)

  # Read back and check the array schema against the tiledb create options
  arr <- tiledb::tiledb_array(uri)
  tsch <- tiledb::schema(arr)

  expect_equal(tiledb::capacity(tsch), 8000)
  expect_equal(tiledb::tile_order(tsch), "COL_MAJOR")
  expect_equal(tiledb::cell_order(tsch), "UNORDERED")

  offsets_filters <- tiledb::filter_list(tsch)$offsets
  expect_equal(tiledb::nfilters(offsets_filters), 1)
  o1 <- offsets_filters[0] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(o1), "RLE")

  validity_filters <- tiledb::filter_list(tsch)$validity
  expect_equal(tiledb::nfilters(validity_filters), 2)
  v1 <- validity_filters[0] # C++ indexing here
  v2 <- validity_filters[1] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(v1), "RLE")
  expect_equal(tiledb::tiledb_filter_type(v2), "NONE")

  dom <- tiledb::domain(tsch)
  expect_equal(tiledb::tiledb_ndim(dom), 2)
  dim0 <- tiledb::dimensions(dom)[[1]]
  expect_equal(tiledb::name(dim0), "soma_dim_0")
  # TODO: As noted above, check this when we are able to.
  # expect_equal(tiledb::tile(dim0), 999)
  dim0_filters <- tiledb::filter_list(dim0)
  expect_equal(tiledb::nfilters(dim0_filters), 3)
  d1 <- dim0_filters[0] # C++ indexing here
  d2 <- dim0_filters[1] # C++ indexing here
  d3 <- dim0_filters[2] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(d1), "RLE")
  expect_equal(tiledb::tiledb_filter_type(d2), "ZSTD")
  expect_equal(tiledb::tiledb_filter_type(d3), "NONE")
  expect_equal(tiledb::tiledb_filter_get_option(d2, "COMPRESSION_LEVEL"), 8)

  dim1 <- tiledb::dimensions(dom)[[2]]
  expect_equal(tiledb::name(dim1), "soma_dim_1")
  # TODO: As noted above, check this when we are able to.
  # expect_equal(tiledb::tile(dim1), 999)
  dim1_filters <- tiledb::filter_list(dim1)
  expect_equal(tiledb::nfilters(dim1_filters), 1)
  d1 <- dim1_filters[0] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(d1), "RLE")

  expect_equal(length(tiledb::attrs(tsch)), 1)
  soma_data_filters <- tiledb::filter_list(tiledb::attrs(tsch)$soma_data)
  expect_equal(tiledb::nfilters(soma_data_filters), 2)

  a1 <- soma_data_filters[0] # C++ indexing here
  a2 <- soma_data_filters[1] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(a1), "BITSHUFFLE")
  expect_equal(tiledb::tiledb_filter_type(a2), "ZSTD")
  expect_equal(tiledb::tiledb_filter_get_option(a2, "COMPRESSION_LEVEL"), 9)

  snda$close()
})

test_that("platform_config defaults", {
  uri <- withr::local_tempdir("soma-sparse-nd-array")

  # Set tiledb create options
  cfg <- PlatformConfig$new()

  # Create the SOMASparseNDArray
  snda <- SOMASparseNDArrayCreate(uri=uri, type=arrow::int32(), shape=c(100,100), platform_config = cfg)

  # Read back and check the array schema against the tiledb create options
  arr <- tiledb::tiledb_array(uri)
  tsch <- tiledb::schema(arr)

  # Here we're snooping on the default dim filter that's used when no other is specified.
  dom <- tiledb::domain(tsch)
  expect_equal(tiledb::tiledb_ndim(dom), 2)

  dim0 <- tiledb::dimensions(dom)[[1]]
  expect_equal(tiledb::name(dim0), "soma_dim_0")
  dim0_filters <- tiledb::filter_list(dim0)
  expect_equal(tiledb::nfilters(dim0_filters), 1)
  d1 <- dim0_filters[0] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(d1), "ZSTD")
  expect_equal(tiledb::tiledb_filter_get_option(d1, "COMPRESSION_LEVEL"), 3)

  dim1 <- tiledb::dimensions(dom)[[2]]
  expect_equal(tiledb::name(dim1), "soma_dim_1")
  dim1_filters <- tiledb::filter_list(dim1)
  expect_equal(tiledb::nfilters(dim1_filters), 1)
  d1 <- dim1_filters[0] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(d1), "ZSTD")
  expect_equal(tiledb::tiledb_filter_get_option(d1, "COMPRESSION_LEVEL"), 3)

  snda$close()
})

test_that("SOMASparseNDArray timestamped ops", {
  uri <- withr::local_tempdir("soma-sparse-nd-array-timestamps")

  # t=10: create 2x2 array and write 1 into top-left entry
  t10 <- Sys.time()
  snda <- SOMASparseNDArrayCreate(uri=uri, type=arrow::int16(), shape=c(2,2))
  snda$write(Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2)))
  snda$close()
  Sys.sleep(1.0)

  # t=20: write 1 into bottom-right entry
  t20 <- Sys.time()
  snda <- SOMASparseNDArrayOpen(uri=uri, mode="WRITE")
  snda$write(Matrix::sparseMatrix(i = 2, j = 2, x = 1, dims = c(2, 2)))
  snda$close()

  # read with no timestamp args and see both writes
  snda <- SOMASparseNDArrayOpen(uri=uri)
  expect_equal(sum(snda$read()$sparse_matrix()$concat()), 2)
  snda$close()

  # read @ t=15 and see only the first write
  snda <- SOMASparseNDArrayOpen(uri=uri, tiledb_timestamp = t10 + 0.5*as.numeric(t20 - t10))
  expect_equal(sum(snda$read()$sparse_matrix()$concat()), 1)
  snda$close()
})
