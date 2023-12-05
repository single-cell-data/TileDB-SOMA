test_that("SOMASparseNDArray creation", {
  skip_if(!extended_tests())
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

  # Verify the array is still open for write
  expect_equal(ndarray$mode(), "WRITE")
  expect_true(tiledb::tiledb_array_is_open(ndarray$object))
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

test_that("SOMASparseNDArray read_sparse_matrix", {
  skip_if(!extended_tests())
  uri <- withr::local_tempdir("sparse-ndarray-3")
  ndarray <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(10, 10))

  # For this test, write 9x9 data into 10x10 array. Leaving the last row & column
  # empty touches corner cases with setting dims() correctly
  mat <- create_sparse_matrix_with_int_dims(10, 10)
  ndarray$write(mat)
  expect_equal(as.numeric(ndarray$shape()), c(10, 10))
  ndarray$close()

  # read_sparse_matrix
  ndarray <- SOMASparseNDArrayOpen(uri)
  mat2 <- ndarray$read()$sparse_matrix(zero_based = T)$concat()
  expect_true(inherits(mat2, "matrixZeroBasedView"))
  expect_s4_class(mat2$get_one_based_matrix(), "sparseMatrix")
  expect_equal(mat2$dim(), c(10, 10))
  expect_equal(mat2$nrow(), 10)
  expect_equal(mat2$ncol(), 10)
  ## not sure why all.equal(mat, mat2) does not pass
  expect_true(all.equal(as.numeric(mat[1:9, 1:9]), as.numeric(mat2$take(0:8, 0:8)$get_one_based_matrix())))
  expect_equal(sum(mat), sum(mat2$get_one_based_matrix()))

  ndarray <- SOMASparseNDArrayOpen(uri)

  ndarray$close()
})

test_that("SOMASparseNDArray read_sparse_matrix_zero_based", {
  skip_if(!extended_tests())
  uri <- withr::local_tempdir("sparse-ndarray")
  ndarray <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(10, 10))

  # For this test, write 9x9 data into 10x10 array. Leaving the last row & column
  # empty touches corner cases with setting dims() correctly
  mat <- create_sparse_matrix_with_int_dims(9, 9)
  ndarray$write(mat)
  expect_equal(as.numeric(ndarray$shape()), c(10, 10))
  ndarray$close()

  # read_sparse_matrix
  ndarray <- SOMASparseNDArrayOpen(uri)
  mat2 <- ndarray$read()$sparse_matrix(zero_based=T)$concat()
  expect_true(inherits(mat2, "matrixZeroBasedView"))
  expect_s4_class(mat2$get_one_based_matrix(), "sparseMatrix")
  expect_equal(mat2$dim(), c(10, 10))
  expect_equal(mat2$nrow(), 10)
  expect_equal(mat2$ncol(), 10)
  ## not sure why all.equal(mat, mat2) does not pass
  expect_true(all.equal(as.numeric(mat), as.numeric(mat2$take(0:8,0:8)$get_one_based_matrix())))
  expect_equal(sum(mat), sum(mat2$get_one_based_matrix()))

  ndarray <- SOMASparseNDArrayOpen(uri)

  # repeat with iterated reader
  iterator <- ndarray$read()$sparse_matrix(zero_based = T)
  mat2 <- iterator$read_next()
  expect_true(inherits(mat2, "matrixZeroBasedView"))
  expect_s4_class(mat2$get_one_based_matrix(), "sparseMatrix")
  expect_equal(mat2$dim(), c(10, 10))
  expect_equal(mat2$nrow(), 10)
  expect_equal(mat2$ncol(), 10)
  expect_true(all.equal(as.numeric(mat), as.numeric(mat2$take(0:8,0:8)$get_one_based_matrix())))
  expect_equal(sum(mat), sum(mat2$get_one_based_matrix()))
  ndarray$close()
})

test_that("SOMASparseNDArray creation with duplicates", {
  skip_if(!extended_tests())
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
  skip_if(!extended_tests())
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
  skip_if(!extended_tests())
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
  skip_if(!extended_tests())
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

test_that("SOMASparseNDArray compatibility with shape >= 2^31 - 1", {
  skip_if(!extended_tests())
  uri <- create_and_populate_32bit_sparse_nd_array(
    uri = withr::local_tempdir("soma-32bit-sparse-nd-array")
  )

  # Coords for all non-zero entries in the array
  all_coords <- bit64::as.integer64(c(0, 2^31 - 2, 2^31 - 1))
  # Coords within R Matrix limits
  safe_coords <- all_coords[1:2]

  snda <- SOMASparseNDArrayOpen(uri, mode = "READ")

  expect_silent(snda$read())
  expect_silent(snda$read()$tables())

  # Arrow table contains all data
  tbl <- snda$read()$tables()$concat()
  expect_identical(tbl$soma_data$as_vector(), c(1L, 2L, 3L))
  expect_identical(tbl$soma_dim_0$as_vector(), as.integer(all_coords))

  # Warning upon creation of SparseReadIter
  expect_warning(
    snda_reader <- snda$read()$sparse_matrix(),
    "Array's shape exceeds"
  )

  # Error when attempting to create a sparse matrix with coordinates >= 2^31-1
  expect_error(
    snda_reader$concat(),
    "Query contains 0-based coordinates outside"
  )

  # Sparse matrix can be created from coordinates within [0, 2^31 - 1]
  suppressWarnings(
    mat <- snda$read(list(safe_coords, safe_coords))$sparse_matrix()$concat()
  )
  expect_identical(dim(mat), as.integer(c(2^31 - 1, 2^31 - 1)))
  expect_length(mat@i, 2)
})

test_that("SOMASparseNDArray bounding box", {
  uri <- withr::local_tempdir("sparse-ndarray-bbox")
  nrows <- 100L
  ncols <- 500L
  ndarray <- SOMASparseNDArrayCreate(uri, type = arrow::int32(), shape = c(nrows, ncols))

  mat <- create_sparse_matrix_with_int_dims(nrows, ncols)
  ndarray$write(mat)
  ndarray$close()

  ndarray <- SOMASparseNDArrayOpen(uri)
  dnames <- ndarray$dimnames()
  bbox_names <- vector('character', length(dnames) * 2L)
  for (i in seq_along(bbox_names)) {
    type <- c('_upper', '_lower')[(i %% 2) + 1L]
    bbox_names[i] <- paste0(dnames[ceiling(i / 2)], '_domain', type)
  }

  expect_true(all(bbox_names %in% names(tiledb::tiledb_get_all_metadata(ndarray$object))))
  for (i in seq_along(bbox_names)) {
    expect_s3_class(x <- ndarray$get_metadata(bbox_names[i]), 'integer64')
    if (i %% 2) {
      expect_equal(x, bit64::as.integer64(0L))
    } else {
      expect_equal(x, bit64::as.integer64(dim(mat)[ceiling(i / 2)] - 1L))
    }
  }

  expect_type(bbox <- ndarray$used_shape(index1 = TRUE), 'list')
  expect_length(bbox, length(dim(mat)))
  expect_equal(names(bbox), dnames)
  expect_true(all(vapply(bbox, length, integer(1L)) == 2L))
  for (i in seq_along(bbox)) {
    expect_equal(bbox[[i]], bit64::as.integer64(c(1L, dim(mat)[i])))
  }

  expect_type(bbox0 <- ndarray$used_shape(index1 = FALSE), 'list')
  expect_length(bbox0, length(dim(mat)))
  expect_equal(names(bbox0), dnames)
  expect_true(all(vapply(bbox0, length, integer(1L)) == 2L))
  for (i in seq_along(bbox0)) {
    expect_equal(bbox0[[i]], bit64::as.integer64(c(0L, dim(mat)[i] - 1L)))
  }

  expect_s3_class(bboxS <- ndarray$used_shape(simplify = TRUE), 'integer64')
  expect_length(bboxS, length(dim(mat)))
  expect_equal(names(bboxS), dnames)
  for (i in seq_along(bboxS)) {
    # Use [[ to remove name from sliced vector
    expect_equal(bboxS[[i]], bit64::as.integer64(dim(mat)[i] - 1L))
  }
})

test_that("SOMASparseNDArray without bounding box", {
  uri <- withr::local_tempdir("sparse-ndarray-no-bbox")
  nrows <- 100L
  ncols <- 500L
  ndarray <- SOMASparseNDArrayCreate(uri, type = arrow::int32(), shape = c(nrows, ncols))

  ndarray$close()

  ndarray <- SOMASparseNDArrayOpen(uri)
  dnames <- ndarray$dimnames()
  bbox_names <- vector('character', length(dnames) * 2L)
  for (i in seq_along(bbox_names)) {
    type <- c('_upper', '_lower')[(i %% 2) + 1L]
    bbox_names[i] <- paste0(dnames[ceiling(i / 2)], '_domain', type)
  }

  expect_false(all(bbox_names %in% names(tiledb::tiledb_get_all_metadata(ndarray$object))))

  expect_error(ndarray$used_shape())
})

test_that("SOMASparseNDArray with failed bounding box", {
  uri <- withr::local_tempdir("sparse-ndarray-failed-bbox")
  nrows <- 100L
  ncols <- 500L
  ndarray <- SOMASparseNDArrayCreate(uri, type = arrow::int32(), shape = c(nrows, ncols))

  mat <- create_sparse_matrix_with_int_dims(nrows, ncols, repr = "T")
  coo <- data.frame(
    i = bit64::as.integer64(slot(mat, "i")),
    j = bit64::as.integer64(slot(mat, "j")),
    x = slot(mat, "x")
  )
  names(coo) <- c(ndarray$dimnames(), ndarray$attrnames())
  ndarray$.__enclos_env__$private$.write_coo_dataframe(coo)

  ndarray$close()

  ndarray <- SOMASparseNDArrayOpen(uri)
  dnames <- ndarray$dimnames()
  bbox_names <- vector('character', length(dnames) * 2L)
  for (i in seq_along(bbox_names)) {
    type <- c('_upper', '_lower')[(i %% 2) + 1L]
    bbox_names[i] <- paste0(dnames[ceiling(i / 2)], '_domain', type)
  }

  expect_false(all(bbox_names %in% names(tiledb::tiledb_get_all_metadata(ndarray$object))))

  expect_error(ndarray$used_shape())
})

test_that("SOMASparseNDArray bounding box implicitly-stored values", {
  uri <- withr::local_tempdir("sparse-ndarray-bbox-implicit")
  nrows <- 100L
  ncols <- 500L
  ndarray <- SOMASparseNDArrayCreate(uri, type = arrow::int32(), shape = c(nrows, ncols))

  mat <- create_sparse_matrix_with_int_dims(nrows, ncols)
  mat[1, ] <- mat[nrows, ] <- mat[, 1L] <- mat[, ncols] <- 0
  ndarray$write(mat)
  ndarray$close()

  ndarray <- SOMASparseNDArrayOpen(uri)
  dnames <- ndarray$dimnames()
  bbox_names <- vector('character', length(dnames) * 2L)
  for (i in seq_along(bbox_names)) {
    type <- c('_upper', '_lower')[(i %% 2) + 1L]
    bbox_names[i] <- paste0(dnames[ceiling(i / 2)], '_domain', type)
  }

  expect_true(all(bbox_names %in% names(tiledb::tiledb_get_all_metadata(ndarray$object))))
  for (i in seq_along(bbox_names)) {
    expect_s3_class(x <- ndarray$get_metadata(bbox_names[i]), 'integer64')
    if (i %% 2) {
      expect_equal(x, bit64::as.integer64(0L))
    } else {
      expect_equal(x, bit64::as.integer64(dim(mat)[ceiling(i / 2)] - 1L))
    }
  }

  expect_type(bbox <- ndarray$used_shape(index1 = TRUE), 'list')
  expect_length(bbox, length(dim(mat)))
  expect_equal(names(bbox), dnames)
  expect_true(all(vapply(bbox, length, integer(1L)) == 2L))
  for (i in seq_along(bbox)) {
    expect_equal(bbox[[i]], bit64::as.integer64(c(1L, dim(mat)[i])))
  }

  expect_type(bbox0 <- ndarray$used_shape(index1 = FALSE), 'list')
  expect_length(bbox0, length(dim(mat)))
  expect_equal(names(bbox0), dnames)
  expect_true(all(vapply(bbox0, length, integer(1L)) == 2L))
  for (i in seq_along(bbox0)) {
    expect_equal(bbox0[[i]], bit64::as.integer64(c(0L, dim(mat)[i] - 1L)))
  }

  expect_s3_class(bboxS <- ndarray$used_shape(simplify = TRUE), 'integer64')
  expect_length(bboxS, length(dim(mat)))
  expect_equal(names(bboxS), dnames)
  for (i in seq_along(bboxS)) {
    # Use [[ to remove name from sliced vector
    expect_equal(bboxS[[i]], bit64::as.integer64(dim(mat)[i] - 1L))
  }

  ranges <- bit64::integer64(2L)
  for (i in seq_along(ranges)) {
    s <- c('i', 'j')[i]
    ranges[i] <- bit64::as.integer64(max(range(slot(mat, s))))
  }
  expect_equal(ndarray$non_empty_domain(), ranges)
  expect_true(all(ndarray$non_empty_domain() < ndarray$used_shape(simplify = TRUE)))
})

test_that("Bounding box assertions", {
  uri <- withr::local_tempdir("bbox-assertions")
  nrows <- 100L
  ncols <- 500L
  ndarray <- SOMASparseNDArrayCreate(uri, type = arrow::int32(), shape = c(nrows, ncols))
  on.exit(ndarray$close())

  mat <- create_sparse_matrix_with_int_dims(nrows, ncols)

  expect_error(ndarray$write(mat, bbox = TRUE))
  expect_error(ndarray$write(mat, bbox = c(TRUE, TRUE)))
  expect_error(ndarray$write(mat, bbox = c(nrows, ncols) + 0.1))
  expect_error(ndarray$write(mat, bbox = list(nrows, ncols)))
  expect_error(ndarray$write(mat, bbox = list(TRUE, TRUE)))
  expect_error(ndarray$write(mat, bbox = c(a = nrows, b = ncols)))
  expect_error(ndarray$write(mat, bbox = list(c(TRUE, FALSE), c('a', 'b'))))
  expect_error(ndarray$write(mat, bbox = c(nrows, ncols) / 2L))
  expect_error(ndarray$write(mat, bbox = list(c(20L, nrows), c(20L, ncols))))
  expect_error(ndarray$write(mat, bbox = list(c(-20L, nrows), c(-20L, ncols))))
})
