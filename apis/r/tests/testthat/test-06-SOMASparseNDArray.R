test_that("SOMASparseNDArray creation", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "sparse-ndarray")
  ndarray <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(10, 10))

  expect_match(
    get_tiledb_object_type(
      ndarray$uri,
      ndarray$.__enclos_env__$private$.soma_context$handle
    ),
    "ARRAY"
  )
  expect_equal(ndarray$dimnames(), c("soma_dim_0", "soma_dim_1"))

  expect_equal(ndarray$attrnames(), "soma_data")
  expect_equal(ndarray$attributes()$soma_data$type, "INT32")

  mat <- create_sparse_matrix_with_int_dims(10, 10)
  vals <- as.vector(t(as.matrix(mat)))
  vals <- vals[vals != 0] # needed below for comparison
  ndarray$write(mat)

  # Verify the array is still open for write
  expect_equal(ndarray$mode(), "WRITE")
  ndarray$close()

  ndarray <- SOMASparseNDArrayOpen(uri)

  # Array write should fail if array opened in read mode
  expect_error(ndarray$write(mat))

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
    coords = list(soma_dim_0 = 0, soma_dim_1 = 0:2),
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
  expect_true(ndarray$is_sparse())
  expect_false(ndarray$allows_duplicates())

  expect_equal(ndarray$shape(), bit64::as.integer64(c(10, 10)))

  expect_true(ndarray$tiledbsoma_has_upgraded_shape())
  shape <- ndarray$shape()
  maxshape <- ndarray$maxshape()
  expect_equal(length(shape), length(maxshape))
  for (i in 1:length(shape)) {
    expect_true(maxshape[i] >= shape[i])
  }

  ## ndim
  expect_equal(ndarray$ndim(), 2L)

  ## nnz
  expect_equal(ndarray$nnz(), 60L)

  ndarray$close()
})

test_that("SOMASparseNDArray write COO assertions", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "sparse-ndarray-coo")
  shape <- c(10L, 10L)
  ndarray <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = shape)

  expect_s3_class(ndarray, "SOMASparseNDArray")
  expect_equal(ndarray$ndim(), 2L)
  mat <- create_sparse_matrix_with_int_dims(10L, 10L)
  df <- data.frame(
    soma_dim_0 = mat@i,
    soma_dim_1 = mat@j,
    soma_data = as.integer(mat@x)
  )

  ndarray$reopen("WRITE")
  expect_invisible(ndarray$.write_coordinates(df))
  ndarray$close()

  # Test write with Table
  tbl <- arrow::as_arrow_table(df)
  ndarray <- SOMASparseNDArrayCreate(
    tempfile(pattern = "sparse-ndarray-coo-table"),
    type = tbl$soma_data$type,
    shape = shape
  )
  expect_invisible(ndarray$.write_coordinates(tbl))
  ndarray$close()

  # Test write unnamed data frame
  udf <- df
  names(udf) <- NULL
  ndarray <- SOMASparseNDArrayCreate(
    uri = tempfile(pattern = "sparse-ndarray-coo-unnamed"),
    type = arrow::int32(),
    shape = shape
  )
  expect_invisible(ndarray$.write_coordinates(udf))
  ndarray$close()

  # Test argument assertions
  arr <- SOMASparseNDArrayCreate(tempfile(), arrow::int32(), shape = shape)
  on.exit(arr$close(), add = TRUE, after = FALSE)
  expect_error(
    arr$.write_coordinates(mat),
    regexp = "'values' must be a data frame or Arrow Table"
  )
  expect_error(
    arr$.write_coordinates(mtcars),
    regexp = "'values' must have one column for each dimension and the data"
  )

  sdf <- df
  while (identical(names(sdf), c(ndarray$dimnames(), ndarray$attrnames()))) {
    sdf <- sdf[, sample(names(sdf)), drop = FALSE]
  }
  expect_error(
    arr$.write_coordinates(sdf),
    regexp = "'values' must be named with the dimension and attribute labels"
  )

  # Test dimension assertions
  ddf <- df
  ddf$soma_dim_0 <- ddf$soma_dim_0 + 0.1
  expect_error(
    arr$.write_coordinates(ddf),
    regexp = "All dimension columns must be integerish"
  )

  ndf <- df
  ndf$soma_dim_0 <- -ndf$soma_dim_0
  expect_error(
    arr$.write_coordinates(ndf),
    regexp = "Dimension columns cannot contain negative values"
  )

  bdf <- df
  bdf$soma_dim_0 <- bdf$soma_dim_0 * 1000
  expect_error(
    arr$.write_coordinates(bdf),
    regexp = "Dimension columns cannot exceed the shape of the array"
  )

  # Test attribute assertions
  ldf <- df
  ldf$soma_data <- TRUE
  expect_error(
    arr$.write_coordinates(ldf),
    regexp = "The data column must be of type 'integer'"
  )

  fdf <- df
  fdf$soma_data <- fdf$soma_data + 0.1
  expect_error(
    arr$.write_coordinates(fdf),
    regexp = "The data column must be of type 'integer'"
  )
})

test_that("SOMASparseNDArray read_sparse_matrix", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "sparse-ndarray-3")
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
  expect_true(all.equal(
    as.numeric(mat[1:9, 1:9]),
    as.numeric(mat2$take(0:8, 0:8)$get_one_based_matrix())
  ))
  expect_equal(sum(mat), sum(mat2$get_one_based_matrix()))

  ndarray <- SOMASparseNDArrayOpen(uri)

  ndarray$close()
})

test_that("SOMASparseNDArray read_sparse_matrix_zero_based", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "sparse-ndarray")
  ndarray <- SOMASparseNDArrayCreate(uri, arrow::int32(), shape = c(10, 10))

  # For this test, write 9x9 data into 10x10 array. Leaving the last row & column
  # empty touches corner cases with setting dims() correctly
  mat <- create_sparse_matrix_with_int_dims(9, 9)
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
  expect_true(all.equal(
    as.numeric(mat),
    as.numeric(mat2$take(0:8, 0:8)$get_one_based_matrix())
  ))
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
  expect_true(all.equal(
    as.numeric(mat),
    as.numeric(mat2$take(0:8, 0:8)$get_one_based_matrix())
  ))
  expect_equal(sum(mat), sum(mat2$get_one_based_matrix()))
  ndarray$close()
})

test_that("SOMASparseNDArray read coordinates", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "sparse-ndarray")
  nrows <- 100L
  ncols <- 20L

  ndarray <- create_and_populate_sparse_nd_array(
    uri = uri,
    mode = "READ",
    nrows = nrows,
    ncols = ncols,
    seed = 42L
  )
  on.exit(ndarray$close(), add = TRUE, after = FALSE)

  expect_identical(as.integer(ndarray$shape()), c(nrows, ncols))
  expect_s4_class(mat <- ndarray$read()$sparse_matrix()$concat(), "dgTMatrix")
  expect_identical(dim(mat), c(nrows, ncols))

  # Note: slices `:` yield integers, not numerics
  # Note: #L is integer, # on its own is numeric
  cases <- list(
    # Test one dim NULL
    "dim0 null, dim1 slice" = list(soma_dim_0 = NULL, soma_dim_1 = 0:9),
    "dim0 null, dim1 slice" = list(soma_dim_0 = 35:45, soma_dim_1 = NULL),
    "dim0 null, dim1 coords" = list(
      soma_dim_0 = NULL,
      soma_dim_1 = c(0L, 5L, 10L)
    ),
    "dim0 coords, dim1 null" = list(soma_dim_0 = c(72, 83), soma_dim_1 = NULL),
    # Test both dims null
    "dim0 null, dim1 null" = list(soma_dim_0 = NULL, soma_dim_1 = NULL),
    # Test both dims provided
    "dim0 coords, dim1 coords" = list(
      soma_dim_0 = c(72, 83),
      soma_dim_1 = c(0L, 5L, 10L)
    ),
    "dim0 slice, dim1 slice" = list(soma_dim_0 = 35:45, soma_dim_1 = 0:9),
    "dim0 coords, dim1 slice" = list(soma_dim_0 = c(72, 83), soma_dim_1 = 0:9),
    "dim0 slice, dim0 coords" = list(
      soma_dim_0 = 35:45,
      soma_dim_1 = c(0L, 5L, 10L)
    ),
    # Test one dim missing
    "dim0 missing, dim1 slice" = list(soma_dim_1 = 0:9),
    "dim0 missing, dim1 coords" = list(soma_dim_1 = c(0L, 5L, 10L)),
    "dim0 missing, dim1 null" = list(soma_dim_1 = NULL),
    "dim0 slice, dim1 missing" = list(soma_dim_0 = 35:45),
    "dim0 coords, dim1 missing" = list(soma_dim_0 = c(72, 83)),
    "dim0 coords, dim1 null" = list(soma_dim_0 = NULL),
    # Test zero-pull
    "zero-pull" = list(soma_dim_0 = c(0, 3), soma_dim_1 = c(0L, 9L))
  )
  for (i in seq_along(cases)) {
    coords <- cases[[i]]
    label <- names(cases)[i]
    expect_s3_class(tbl <- ndarray$read(coords)$tables()$concat(), "Table")
    ii <- if (is.null(coords$soma_dim_0)) {
      TRUE
    } else {
      mat@i %in% coords$soma_dim_0
    }
    jj <- if (is.null(coords$soma_dim_1)) {
      TRUE
    } else {
      mat@j %in% coords$soma_dim_1
    }
    nr <- ifelse(
      isTRUE(ii) && isTRUE(jj),
      yes = length(mat@x),
      no = sum(ii & jj)
    )
    expect_identical(nrow(tbl), nr, label = label)
  }

  # Test assertions
  list_cases <- list(TRUE, "tomato", 1L, 1.1, bit64::as.integer64(1L), list())
  for (coords in list_cases) {
    expect_error(ndarray$read(coords), regexp = "'coords' must be a list")
  }

  intgerish_cases <- list(
    list(TRUE),
    list("tomato"),
    list(1.1),
    list(NA_integer_),
    list(NA_real_),
    list(bit64::NA_integer64_),
    list(Inf),
    list(-4),
    list(factor(letters[1:10])),
    list(matrix(1:10, ncol = 1:10)),
    list(array(1:10))
  )
  for (coords in intgerish_cases) {
    expect_error(
      ndarray$read(coords),
      regexp = "'coords' must be a list integerish vectors"
    )
  }

  names_cases <- list(
    list(1:3, 1:5, 1:10),
    list(tomato = 1:10),
    list(soma_dim_0 = 1:10, tomato = 1:10)
  )
  for (coords in names_cases) {
    expect_error(
      ndarray$read(coords),
      regexp = "'coords' if unnamed must have length of dim names, else if named names must match dim names"
    )
  }
})

test_that("SOMASparseNDArray creation with duplicates", {
  skip_if(!extended_tests())

  set.seed(42)
  D <- data.frame(
    soma_dim_0 = sample(seq.int(0L, 99L), size = 10L),
    soma_dim_1 = sample(seq.int(0L, 99L), size = 10L),
    soma_data = rnorm(n = 10L)
  )
  params <- list(
    c(allow_dups = FALSE, do_dup = FALSE),
    c(allow_dups = TRUE, do_dup = FALSE),
    c(allow_dups = TRUE, do_dup = TRUE)
  )
  for (i in seq_along(along.with = params)) {
    allow <- params[[i]][['allow_dups']]
    dups <- params[[i]][['do_dup']]
    cfg <- if (allow) {
      PlatformConfig$new()$set(
        platform = "tiledb",
        param = "create",
        key = "allows_duplicates",
        value = TRUE
      )
    } else {
      NULL
    }
    arr <- SOMASparseNDArrayCreate(
      uri = tempfile(
        pattern = sprintf(
          fmt = "sparse-ndarray-%s-%s",
          ifelse(allow, yes = "dupAllowed", no = "dupDisallowed"),
          ifelse(dups, yes = "duplicated", no = "nodups")
        )
      ),
      type = arrow::float64(),
      shape = c(100L, 100L),
      platform_config = cfg
    )
    arr$.write_coordinates(if (dups) rbind(D, D) else D)
    arr$reopen("READ")
    expect_equal(arr$nnz(), if (dups) nrow(D) * 2L else nrow(D))
    arr$close()
    unlink(arr$uri, recursive = TRUE, force = TRUE)
  }
})

test_that("platform_config is respected", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-sparse-nd-array")

  # Set tiledb create options
  cfg <- PlatformConfig$new()
  cfg$set("tiledb", "create", "sparse_nd_array_dim_zstd_level", 9)
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

  # Create the SOMASparseNDArray
  snda <- SOMASparseNDArrayCreate(
    uri = uri,
    type = arrow::int32(),
    shape = c(100, 100),
    platform_config = cfg
  )

  # Read back and check the array schema against the tiledb create options
  expect_equal(
    c_capacity(snda$uri, snda$.__enclos_env__$private$.soma_context$handle),
    8000L
  )
  expect_equal(
    c_tile_order(snda$uri, snda$.__enclos_env__$private$.soma_context$handle),
    "COL_MAJOR"
  )
  expect_equal(
    c_cell_order(snda$uri, snda$.__enclos_env__$private$.soma_context$handle),
    "ROW_MAJOR"
  )

  expect_length(
    coord_filters <- c_schema_filters(
      snda$uri,
      snda$.__enclos_env__$private$.soma_context$handle
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
    domain <- c_domain(snda$uri, snda$.__enclos_env__$private$.soma_context$handle),
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
  # expect_equal(tiledb::tile(dim), 999)
  dim0 <- domain$soma_dim_0
  expect_length(dim0$filters, n = 3L)
  expect_equal(dim0$filters[[1L]]$filter_type, "RLE")
  expect_equal(dim0$filters[[2L]]$filter_type, "ZSTD")
  expect_equal(dim0$filters[[2L]]$compression_level, 8L)
  expect_equal(dim0$filters[[3L]]$filter_type, "NOOP")

  dim1 <- domain$soma_dim_1
  expect_length(dim1$filters, n = 1L)
  expect_equal(dim1$filters[[1L]]$filter_type, "RLE")

  expect_length(attrs <- snda$attributes(), n = 1L)
  expect_length(attrs$soma_data$filter_list, n = 2L)
  expect_equal(attrs$soma_data$filter_list[[1L]]$filter_type, "BITSHUFFLE")
  expect_equal(attrs$soma_data$filter_list[[2L]]$filter_type, "ZSTD")
  expect_equal(attrs$soma_data$filter_list[[2L]]$compression_level, 9L)

  snda$close()
})

test_that("platform_config defaults", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-sparse-nd-array")

  # Set tiledb create options
  cfg <- PlatformConfig$new()

  # Create the SOMASparseNDArray
  snda <- SOMASparseNDArrayCreate(
    uri = uri,
    type = arrow::int32(),
    shape = c(100, 100),
    platform_config = cfg
  )

  # Here we're snooping on the default dim filter that's used when no other is specified.
  expect_length(
    domain <- c_domain(snda$uri, snda$.__enclos_env__$private$.soma_context$handle),
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

  snda$close()
})

test_that("SOMASparseNDArray timestamped ops", {
  skip_if(!extended_tests())
  uri <- tempfile(pattern = "soma-sparse-nd-array-timestamps")

  # t=10: create 2x2 array and write 1 into top-left entry
  t10 <- as.POSIXct(10, tz = "UTC", origin = "1970-01-01")
  snda <- SOMASparseNDArrayCreate(
    uri = uri,
    type = arrow::int16(),
    shape = c(2, 2),
    tiledb_timestamp = t10
  )
  snda$write(Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(2, 2)))
  snda$close()

  # t=20: write 1 into bottom-right entry
  t20 <- as.POSIXct(20, tz = "UTC", origin = "1970-01-01")
  snda <- SOMASparseNDArrayOpen(
    uri = uri,
    mode = "WRITE",
    tiledb_timestamp = t20
  )
  snda$write(Matrix::sparseMatrix(i = 2, j = 2, x = 1, dims = c(2, 2)))
  snda$close()

  # read with no timestamp args and see both writes
  snda <- SOMASparseNDArrayOpen(uri = uri)
  expect_equal(sum(snda$read()$sparse_matrix()$concat()), 2)
  snda$close()

  # read @ t=15 and see only the first write
  snda <- SOMASparseNDArrayOpen(
    uri = uri,
    tiledb_timestamp = t10 + 0.5 * as.numeric(t20 - t10)
  )
  expect_equal(sum(snda$read()$sparse_matrix()$concat()), 1)
  snda$close()
})

test_that("SOMASparseNDArray compatibility with shape >= 2^31 - 1", {
  skip_if(!extended_tests())
  uri <- create_and_populate_32bit_sparse_nd_array(
    uri = tempfile(pattern = "soma-32bit-sparse-nd-array")
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
