test_that("SOMADenseNDArray creation", {
  skip_if(!extended_tests())
  uri <- withr::local_tempdir("dense-ndarray")

  ndarray <- SOMADenseNDArrayCreate(uri, arrow::int32(), shape = c(10, 5))

  expect_equal(tiledb::tiledb_object_type(uri), "ARRAY")
  expect_equal(ndarray$dimnames(), c("soma_dim_0", "soma_dim_1"))
  expect_equal(ndarray$attrnames(), "soma_data")
  expect_equal(tiledb::datatype(ndarray$attributes()$soma_data), "INT32")

  mat <- create_dense_matrix_with_int_dims(10, 5)
  ndarray$write(mat)

  # Verify the array is still open for write
  expect_equal(ndarray$mode(), "WRITE")
  expect_true(tiledb::tiledb_array_is_open(ndarray$object))
  ndarray$close()

  # Read result in column-major order to match R matrix layout
  ndarray <- SOMADenseNDArrayOpen(uri)
  tbl <- ndarray$read_arrow_table(result_order = "COL_MAJOR")
  expect_true(is_arrow_table(tbl))
  expect_equal(tbl$ColumnNames(), c("soma_data"))

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

  ndarray$close()
})

test_that("platform_config is respected", {
  skip_if(!extended_tests())
  uri <- withr::local_tempdir("soma-dense-nd-array")

  # Set tiledb create options
  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dense_nd_array_dim_zstd_level', 9)
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

  # Create the SOMADenseNDArray
  dnda <- SOMADenseNDArrayCreate(uri=uri, type=arrow::int32(), shape=c(100,100), platform_config = cfg)

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

  dnda$close()
})

test_that("platform_config defaults", {
  skip_if(!extended_tests())
  uri <- withr::local_tempdir("soma-dense-nd-array")

  # Set tiledb create options
  cfg <- PlatformConfig$new()

  # Create the SOMADenseNDArray
  dnda <- SOMADenseNDArrayCreate(uri=uri, type=arrow::int32(), shape=c(100,100), platform_config = cfg)

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

  dnda$close()
})

test_that("SOMADenseNDArray timestamped ops", {
  skip_if(!extended_tests())
  uri <- withr::local_tempdir("soma-dense-nd-array-timestamps")

  t10 <- Sys.time()
  dnda <- SOMADenseNDArrayCreate(uri=uri, type=arrow::int16(), shape=c(2,2))
  M1 <- matrix(rep(1, 4), 2, 2)
  dnda$write(M1)
  dnda$close()

  dnda <- SOMADenseNDArrayOpen(uri=uri)
  expect_equal(dnda$read_dense_matrix(), M1)
  dnda$close()
  Sys.sleep(1.0)

  t20 <- Sys.time()
  dnda <- SOMADenseNDArrayOpen(uri=uri, mode="WRITE")
  M2 <- matrix(rep(1, 4), 2, 2)
  dnda$write(M2)
  dnda$close()

  dnda <- SOMADenseNDArrayOpen(uri=uri)
  expect_equal(dnda$read_dense_matrix(), M2)
  dnda$close()

  dnda <- SOMADenseNDArrayOpen(uri=uri, tiledb_timestamp = t10 + 0.5*as.numeric(t20 - t10))
  expect_equal(dnda$read_dense_matrix(), M1)   # read between t10 and t20 sees only first write
  dnda$close()
})
