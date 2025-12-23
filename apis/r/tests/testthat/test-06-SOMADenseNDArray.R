test_that("SOMADenseNDArray creation", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "dense-ndarray")

  ndarray <- SOMADenseNDArrayCreate(uri, arrow::int32(), shape = c(10, 5))

  # The array is open for write on create. Nonetheless, close and
  # reopen to ensure that state needed for a write is available.
  ndarray$close()
  ndarray <- SOMADenseNDArrayOpen(uri, "WRITE")

  expect_match(
    get_tiledb_object_type(
      ndarray$uri,
      ndarray$.__enclos_env__$private$.context$handle
    ),
    "ARRAY"
  )
  expect_equal(ndarray$dimnames(), c("soma_dim_0", "soma_dim_1"))
  expect_equal(ndarray$attrnames(), "soma_data")
  expect_equal(ndarray$attributes()$soma_data$type, "INT32")

  mat <- create_dense_matrix_with_int_dims(10, 5)
  ndarray$write(mat)

  # Verify the array is still open for write
  expect_equal(ndarray$mode(), "WRITE")
  ndarray$close()

  # Read result in column-major order to match R matrix layout
  ndarray <- SOMADenseNDArrayOpen(uri)

  # Array write should fail if array opened in read mode
  expect_error(ndarray$write(mat))

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
    coords = list(soma_dim_0 = 0:3, soma_dim_1 = 0:2),
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
    regexp = "'coords' must be a list integerish vectors"
  )

  # Validate TileDB array schema
  expect_false(ndarray$is_sparse())

  ## shape
  expect_equal(ndarray$shape(), bit64::as.integer64(c(10, 5)))

  ## ndim
  expect_equal(ndarray$ndim(), 2L)

  ndarray$close()
})

test_that("platform_config is respected", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-dense-nd-array")

  # Set tiledb create options
  cfg <- PlatformConfig$new()
  cfg$set("tiledb", "create", "dense_nd_array_dim_zstd_level", 9)
  cfg$set("tiledb", "create", "capacity", 8000)
  cfg$set("tiledb", "create", "tile_order", "COL_MAJOR")
  cfg$set("tiledb", "create", "cell_order", "ROW_MAJOR")
  cfg$set("tiledb", "create", "offsets_filters", list("RLE"))
  cfg$set("tiledb", "create", "validity_filters", list("RLE", "NONE"))
  cfg$set(
    "tiledb",
    "create",
    "dims",
    list(
      soma_dim_0 = list(
        filters = list(
          "RLE",
          list(name = "ZSTD", COMPRESSION_LEVEL = 8),
          "NONE"
        )
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
    )
  )
  cfg$set(
    "tiledb",
    "create",
    "attrs",
    list(
      soma_data = list(
        filters = list("BITSHUFFLE", list(name = "ZSTD", COMPRESSION_LEVEL = 9))
      )
    )
  )

  # Create the SOMADenseNDArray
  dnda <- SOMADenseNDArrayCreate(
    uri = uri,
    type = arrow::int32(),
    shape = c(100, 100),
    platform_config = cfg
  )

  # Read back and check the array schema against the tiledb create options
  expect_equal(
    c_capacity(dnda$uri, dnda$.__enclos_env__$private$.context$handle),
    8000L
  )
  expect_equal(
    c_tile_order(dnda$uri, dnda$.__enclos_env__$private$.context$handle),
    "COL_MAJOR"
  )
  expect_equal(
    c_cell_order(dnda$uri, dnda$.__enclos_env__$private$.context$handle),
    "ROW_MAJOR"
  )

  expect_length(
    coord_filters <- c_schema_filters(
      dnda$uri,
      dnda$.__enclos_env__$private$.context$handle
    ),
    n = 3L
  )
  expect_named(coord_filters, c("coords", "offsets", "validity"))

  expect_length(coord_filters$offsets, n = 1L)
  expect_equal(coord_filters$offsets[[1L]]$filter_type, "RLE")

  expect_length(coord_filters$validity, n = 2L)
  expect_equal(coord_filters$validity[[1L]]$filter_type, "RLE")
  expect_equal(coord_filters$validity[[2L]]$filter_type, "NOOP")

  expect_length(
    domain <- c_domain(dnda$uri, dnda$.__enclos_env__$private$.context$handle),
    n = 2L
  )
  expect_named(domain, dims <- sprintf("soma_dim_%i", 0:1))
  expect_equal(
    vapply(
      domain,
      FUN = '[[',
      FUN.VALUE = character(1L),
      "name",
      USE.NAMES = FALSE
    ),
    dims
  )
  # TODO: As noted above, check this when we are able to.
  # expect_equal(tiledb::tile(dim0), 999)
  dim0 <- domain$soma_dim_0
  expect_length(dim0$filters, n = 3L)
  expect_equal(dim0$filters[[1L]]$filter_type, "RLE")
  expect_equal(dim0$filters[[2L]]$filter_type, "ZSTD")
  expect_equal(dim0$filters[[2L]]$compression_level, 8L)
  expect_equal(dim0$filters[[3L]]$filter_type, "NOOP")

  dim1 <- domain$soma_dim_1
  expect_length(dim1$filters, n = 1L)
  expect_equal(dim1$filters[[1L]]$filter_type, "RLE")

  expect_length(attrs <- dnda$attributes(), n = 1L)
  expect_length(attrs$soma_data$filter_list, n = 2L)
  expect_equal(attrs$soma_data$filter_list[[1L]]$filter_type, "BITSHUFFLE")
  expect_equal(attrs$soma_data$filter_list[[2L]]$filter_type, "ZSTD")
  expect_equal(attrs$soma_data$filter_list[[2L]]$compression_level, 9L)

  dnda$close()
})

test_that("platform_config defaults", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-dense-nd-array")

  # Set tiledb create options
  cfg <- PlatformConfig$new()

  # Create the SOMADenseNDArray
  dnda <- SOMADenseNDArrayCreate(
    uri = uri,
    type = arrow::int32(),
    shape = c(100, 100),
    platform_config = cfg
  )

  # Here we're snooping on the default dim filter that's used when no other is specified.
  expect_length(
    domain <- c_domain(dnda$uri, dnda$.__enclos_env__$private$.context$handle),
    n = 2L
  )
  expect_named(domain, dims <- sprintf("soma_dim_%i", 0:1))
  expect_equal(
    vapply(
      domain,
      FUN = '[[',
      FUN.VALUE = character(1L),
      "name",
      USE.NAMES = FALSE
    ),
    dims
  )
  dim0 <- domain$soma_dim_0
  expect_length(dim0$filters, n = 1L)
  expect_equal(dim0$filters[[1L]]$filter_type, "ZSTD")
  expect_equal(dim0$filters[[1L]]$compression_level, 3L)

  dim1 <- domain$soma_dim_1
  expect_length(dim1$filters, n = 1L)
  expect_equal(dim1$filters[[1L]]$filter_type, "ZSTD")
  expect_equal(dim1$filters[[1L]]$compression_level, 3L)

  dnda$close()
})

test_that("SOMADenseNDArray timestamped ops", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-dense-nd-array-timestamps")

  t10 <- Sys.time()
  dnda <- SOMADenseNDArrayCreate(
    uri = uri,
    type = arrow::int16(),
    shape = c(2, 2)
  )
  M1 <- matrix(rep(1, 4), 2, 2)
  dnda$write(M1)
  dnda$close()

  dnda <- SOMADenseNDArrayOpen(uri = uri)
  expect_equal(dnda$read_dense_matrix(), M1)
  dnda$close()
  Sys.sleep(1.0)

  t20 <- Sys.time()
  dnda <- SOMADenseNDArrayOpen(uri = uri, mode = "WRITE")

  M2 <- matrix(rep(1, 4), 2, 2)
  dnda$write(M2)
  dnda$close()

  dnda <- SOMADenseNDArrayOpen(uri = uri)
  expect_equal(dnda$read_dense_matrix(), M2)
  dnda$close()

  dnda <- SOMADenseNDArrayOpen(
    uri = uri,
    tiledb_timestamp = t10 + 0.5 * as.numeric(t20 - t10)
  )
  expect_equal(dnda$read_dense_matrix(), M1) # read between t10 and t20 sees only first write
  dnda$close()
})

test_that("`SOMADenseNDArray$set_data_type()` deprecations", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-dense-nd-array-timestamps")

  dnda <- SOMADenseNDArrayCreate(
    uri = uri,
    type = arrow::int16(),
    shape = c(2, 2)
  )
  M1 <- matrix(rep(1, 4), 2, 2)
  dnda$write(M1)
  dnda$close()

  dnda <- SOMADenseNDArrayOpen(uri = uri)

  with_mocked_bindings(
    .tiledbsoma_deprecation_version = function() "2.1.0",
    .deprecation_stage = function(when) "deprecate",
    {
      lifecycle::expect_deprecated(dnda$set_data_type(arrow::int16()))
    }
  )

  if (utils::packageVersion("tiledbsoma") >= "2.1.0") {
    lifecycle::expect_deprecated(dnda$set_data_type(arrow::int16()))
  }
})
