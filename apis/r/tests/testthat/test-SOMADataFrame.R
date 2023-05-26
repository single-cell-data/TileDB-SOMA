test_that("Basic mechanics", {
  uri <- withr::local_tempdir("soma-dataframe")
  asch <- create_arrow_schema()

  expect_error(
    SOMADataFrameCreate(uri, asch, index_column_names = "qux"),
    "The following field does not exist: qux"
  )

  sdf <- SOMADataFrameCreate(uri, asch, index_column_names = "foo")
  expect_true(sdf$exists())
  expect_true(dir.exists(uri))

  # check for missing columns
  expect_error(
    sdf$write(arrow::arrow_table(foo = 1L:10L)),
    "All schema fields must be present in 'values'"
  )
  # check for extra columns
  expect_error(
    sdf$write(arrow::arrow_table(qux = 1L:10L)),
    "All columns in 'values' must be defined in the schema"
  )

  tbl0 <- arrow::arrow_table(foo = 1L:36L,
                             soma_joinid = 1L:36L,
                             bar = 1.1:36.1,
                             baz = c("á", "ą", "ã", "à", "å", "ä", "æ", "ç", "ć", "Ç", "í",
                                     "ë", "é", "è", "ê", "ł", "Ł", "ñ", "ń", "ó", "ô", "ò",
                                     "ö", "ø", "Ø", "ř", "š", "ś", "ş", "Š", "ú", "ü", "ý",
                                     "ź", "Ž", "Ż"),
                             schema = asch)

  sdf$write(tbl0)
  sdf$close()

  # Read back the data (ignore attributes)
  sdf <- SOMADataFrameOpen(uri)
  expect_match(sdf$soma_type, "SOMADataFrame")

  expect_equivalent(
    tiledb::tiledb_array(sdf$uri, return_as = "asis")[],
    as.list(tbl0),
    ignore_attr = TRUE
  )

  # Read result should recreate the original Table
  tbl1 <- sdf$read()$concat()
  expect_true(tbl1$Equals(tbl0))
  sdf$close()

  # Same as above but now for RecordBatch
  sdf <- SOMADataFrameOpen(uri, mode = "WRITE")
  rb0 <- arrow::record_batch(foo = 1L:36L,
                             soma_joinid = 1L:36L,
                             bar = 1.1:36.1,
                             baz = c("á", "ą", "ã", "à", "å", "ä", "æ", "ç", "ć", "Ç", "í",
                                     "ë", "é", "è", "ê", "ł", "Ł", "ñ", "ń", "ó", "ô", "ò",
                                     "ö", "ø", "Ø", "ř", "š", "ś", "ş", "Š", "ú", "ü", "ý",
                                     "ź", "Ž", "Ż"),
                             schema = asch)

  sdf$write(rb0)
  sdf$close()

  # Read back the data (ignore attributes)
  sdf <- SOMADataFrameOpen(uri)
  expect_equivalent(
    tiledb::tiledb_array(sdf$uri, return_as = "asis")[],
    as.list(rb0),
    ignore_attr = TRUE
  )

  # Read result should recreate the original RecordBatch (when seen as a tibble)
  rb1 <- arrow::as_record_batch(sdf$read()$concat())
  expect_equivalent(dplyr::collect(rb0), dplyr::collect(rb1))

  # Slicing by foo
  tbl1 <- sdf$read(coords = list(foo = 1L:2L))$concat()
  expect_true(tbl1$Equals(tbl0$Slice(offset = 0, length = 2)))

  # Slicing unnamed also work
  tbl1 <- sdf$read(coords = 1L:2L)$concat()
  expect_true(tbl1$Equals(tbl0$Slice(offset = 0, length = 2)))

  # Subselecting columns
  expect_error(
    sdf$read(column_names = "foobar"),
    "'column_names' must only contain valid dimension or attribute columns"
  )

  tbl1 <- sdf$read(column_names = "bar")$concat()
  expect_true(tbl1$Equals(tbl0$SelectColumns(2L)))

  # Attribute filters
  tbl1 <- sdf$read(value_filter = "bar < 5")$concat()
  expect_true(tbl1$Equals(tbl0$Filter(tbl0$bar < 5)))

  # Validate TileDB array schema
  arr <- tiledb::tiledb_array(uri)
  sch <- tiledb::schema(arr)
  expect_true(tiledb::is.sparse(sch))
  expect_false(tiledb::allows_dups(sch))
  sdf$close()
})

test_that("Basic mechanics with default index_column_names", {
  uri <- withr::local_tempdir("soma-dataframe-soma-joinid")
  asch <- create_arrow_schema(foo_first=FALSE)

  sdf <- SOMADataFrame$new(uri, internal_use_only = "allowed_use")
  expect_error(
    sdf$create(asch, index_column_names = "qux", internal_use_only = "allowed_use"),
    "The following field does not exist: qux"
  )

  sdf$create(asch, internal_use_only = "allowed_use")
  expect_true(sdf$exists())
  expect_true(dir.exists(uri))
  expect_match(sdf$soma_type, "SOMADataFrame")

  # check for missing columns
  expect_error(
    sdf$write(arrow::arrow_table(foo = 1L:10L)),
    "All schema fields must be present in 'values'"
  )
  # check for extra columns
  expect_error(
    sdf$write(arrow::arrow_table(qux = 1L:10L)),
    "All columns in 'values' must be defined in the schema"
  )

  tbl0 <- arrow::arrow_table( soma_joinid = 1L:36L,
                             foo = 1L:36L,
                             bar = 1.1:36.1,
                             baz = c("á", "ą", "ã", "à", "å", "ä", "æ", "ç", "ć", "Ç", "í",
                                     "ë", "é", "è", "ê", "ł", "Ł", "ñ", "ń", "ó", "ô", "ò",
                                     "ö", "ø", "Ø", "ř", "š", "ś", "ş", "Š", "ú", "ü", "ý",
                                     "ź", "Ž", "Ż"),
                             schema = asch)

  sdf$write(tbl0)

  # read back the data (ignore attributes)
  expect_equivalent(
    tiledb::tiledb_array(sdf$uri, return_as = "asis")[],
    as.list(tbl0),
    ignore_attr = TRUE
  )
})

test_that("creation with all supported dimension data types", {

  sch <- arrow::schema(
    arrow::field("int8", arrow::int8(), nullable = FALSE),
    arrow::field("int16", arrow::int16(), nullable = FALSE),
    arrow::field("double", arrow::float64(), nullable = FALSE),
    arrow::field("int", arrow::int32(), nullable = FALSE),
    arrow::field("int64", arrow::int64(), nullable = FALSE),
    arrow::field("string", arrow::utf8(), nullable = FALSE)
  )

  tbl0 <- arrow::arrow_table(
    int8 = 1L:36L,
    int16 = 1:36L,
    double = 1.1:36.1,
    int = 1L:36L,
    int64 = bit64::as.integer64(1L:36L),
    string = c("á", "ą", "ã", "à", "å", "ä", "æ", "ç", "ć", "Ç", "í",
               "ë", "é", "è", "ê", "ł", "Ł", "ñ", "ń", "ó", "ô", "ò",
               "ö", "ø", "Ø", "ř", "š", "ś", "ş", "Š", "ú", "ü", "ý",
               "ź", "Ž", "Ż"),
    schema = sch
  )

  for (dtype in tbl0$ColumnNames()) {
    uri <- withr::local_tempdir(paste0("soma-dataframe-", dtype))
    expect_silent(
      sdf <- SOMADataFrameCreate(uri, tbl0$schema, index_column_names = dtype)
    )
    expect_true(sdf$exists())
    sdf$close()
  }
})

test_that("int64 values are stored correctly", {
  uri <- withr::local_tempdir("soma-dataframe")
  asch <- arrow::schema(
    arrow::field("foo", arrow::int32(), nullable = FALSE),
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
  )

  sdf <- SOMADataFrameCreate(uri, asch, index_column_names = "foo")
  tbl0 <- arrow::arrow_table(foo = 1L:10L, soma_joinid = 1L:10L, schema = asch)

  orig_downcast_value <- getOption("arrow.int64_downcast")

  sdf$write(tbl0)
  sdf$close()

  sdf <- SOMADataFrameOpen(uri)
  tbl1 <- sdf$read()$concat()
  expect_true(tbl1$Equals(tbl0))

  # verify int64_downcast option was restored
  expect_equal(getOption("arrow.int64_downcast"), orig_downcast_value)
  sdf$close()
})

test_that("SOMADataFrame read", {
    uri <- extract_dataset("soma-dataframe-pbmc3k-processed-obs")

    sdf <- SOMADataFrameOpen(uri)
    z <- sdf$read()$concat()
    expect_equal(z$num_rows, 2638L)
    expect_equal(z$num_columns, 6L)
    sdf$close()

    columns <- c("n_counts", "n_genes", "louvain")
    sdf <- SOMADataFrameOpen(uri)
    z <- sdf$read(column_names=columns)$concat()
    expect_equal(z$num_columns, 3L)
    expect_equal(z$ColumnNames(), columns)
    sdf$close()

    columns <- c("n_counts", "does_not_exist")
    sdf <- SOMADataFrameOpen(uri)
    expect_error(sdf$read(column_names=columns))
    sdf$close()

    coords <- bit64::as.integer64(seq(100, 109))
    sdf <- SOMADataFrameOpen(uri)
    z <- sdf$read(coords = list(soma_joinid=coords))$concat()
    expect_equal(z$num_rows, 10L)
    sdf$close()
})

test_that("soma_ prefix is reserved", {
  uri <- withr::local_tempdir("soma-dataframe")
  asch <- create_arrow_schema()

  # Add a soma_joinid column with the wrong type
  asch <- asch$AddField(
    i = 1,
    field = arrow::field("soma_foo", arrow::int32(), nullable = FALSE)
  )

  expect_error(
    SOMADataFrameCreate(uri, asch, index_column_names = "foo"),
    "Column names must not start with reserved prefix 'soma_'"
  )
})

test_that("soma_joinid is added on creation", {
  uri <- withr::local_tempdir("soma-dataframe")
  asch <- create_arrow_schema()
  asch <- asch$RemoveField(match("soma_joinid", asch$names) - 1)

  sdf <- SOMADataFrameCreate(uri, asch, index_column_names = "foo")

  expect_true("soma_joinid" %in% sdf$attrnames())
  expect_equal(tiledb::datatype(sdf$attributes()$soma_joinid), "INT64")
  sdf$close()
})

test_that("soma_joinid validations", {
  uri <- withr::local_tempdir("soma-dataframe")
  asch <- create_arrow_schema()

  # Add a soma_joinid column with the wrong type
  asch <- asch$RemoveField(match("soma_joinid", asch$names) - 1)
  asch <- asch$AddField(
    i = 1,
    field = arrow::field("soma_joinid", arrow::int32(), nullable = FALSE)
  )

  expect_error(
    SOMADataFrameCreate(uri, asch, index_column_names = "foo"),
    "soma_joinid field must be of type Arrow int64"
  )
})

test_that("platform_config is respected", {
  uri <- withr::local_tempdir("soma-dataframe")

  # Set Arrow schema
  asch <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("i32", arrow::int32(), nullable = FALSE),
    arrow::field("f64", arrow::float64(), nullable = FALSE),
    arrow::field("utf8", arrow::large_utf8(), nullable = FALSE)
  )

  # Set tiledb create options
  cfg <- PlatformConfig$new()
  cfg$set('tiledb', 'create', 'dataframe_dim_zstd_level', 8)
  cfg$set('tiledb', 'create', 'sparse_nd_array_dim_zstd_level', 9)
  cfg$set('tiledb', 'create', 'capacity', 8000)
  cfg$set('tiledb', 'create', 'tile_order', 'COL_MAJOR')
  cfg$set('tiledb', 'create', 'cell_order', 'UNORDERED')
  cfg$set('tiledb', 'create', 'offsets_filters', list("RLE"))
  cfg$set('tiledb', 'create', 'validity_filters', list("RLE", "NONE"))
  cfg$set('tiledb', 'create', 'dims', list(
    soma_joinid = list(
      filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=8), "NONE")
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
    i32 = list(
      filters = list("RLE", list(name="ZSTD", COMPRESSION_LEVEL=9))
    ),
    f64 = list(
      filters = list()
    )
  ))

  # Create the SOMADataFrame
  sdf <- SOMADataFrameCreate(uri=uri, schema=asch, index_column_names=c("soma_joinid"), platform_config = cfg)

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
  expect_equal(tiledb::tiledb_ndim(dom), 1)
  dim <- tiledb::dimensions(dom)[[1]]
  expect_equal(tiledb::name(dim), "soma_joinid")
  # TODO: As noted above, check this when we are able to.
  # expect_equal(tiledb::tile(dim), 999)
  dim_filters <- tiledb::filter_list(dim)
  expect_equal(tiledb::nfilters(dim_filters), 3)
  d1 <- dim_filters[0] # C++ indexing here
  d2 <- dim_filters[1] # C++ indexing here
  d3 <- dim_filters[2] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(d1), "RLE")
  expect_equal(tiledb::tiledb_filter_type(d2), "ZSTD")
  expect_equal(tiledb::tiledb_filter_type(d3), "NONE")
  expect_equal(tiledb::tiledb_filter_get_option(d2, "COMPRESSION_LEVEL"), 8)

  expect_equal(length(tiledb::attrs(tsch)), 3)
  i32_filters <- tiledb::filter_list(tiledb::attrs(tsch)$i32)
  f64_filters <- tiledb::filter_list(tiledb::attrs(tsch)$f64)
  expect_equal(tiledb::nfilters(i32_filters), 2)
  expect_equal(tiledb::nfilters(f64_filters), 0)

  i1 <- i32_filters[0] # C++ indexing here
  i2 <- i32_filters[1] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(i1), "RLE")
  expect_equal(tiledb::tiledb_filter_type(i2), "ZSTD")
  expect_equal(tiledb::tiledb_filter_get_option(i2, "COMPRESSION_LEVEL"), 9)

  sdf$close()
})

test_that("platform_config defaults", {
  uri <- withr::local_tempdir("soma-dataframe")

  # Set Arrow schema
  asch <- arrow::schema(
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
    arrow::field("i32", arrow::int32(), nullable = FALSE),
    arrow::field("f64", arrow::float64(), nullable = FALSE),
    arrow::field("utf8", arrow::large_utf8(), nullable = FALSE)
  )

  # Set tiledb create options
  cfg <- PlatformConfig$new()

  # Create the SOMADataFrame
  sdf <- SOMADataFrameCreate(uri=uri, schema=asch, index_column_names=c("soma_joinid"), platform_config = cfg)

  # Read back and check the array schema against the tiledb create options
  arr <- tiledb::tiledb_array(uri)
  tsch <- tiledb::schema(arr)

  # Here we're snooping on the default dim filter that's used when no other is specified.
  dom <- tiledb::domain(tsch)
  expect_equal(tiledb::tiledb_ndim(dom), 1)
  dim <- tiledb::dimensions(dom)[[1]]
  expect_equal(tiledb::name(dim), "soma_joinid")
  dim_filters <- tiledb::filter_list(dim)
  expect_equal(tiledb::nfilters(dim_filters), 1)
  d1 <- dim_filters[0] # C++ indexing here
  expect_equal(tiledb::tiledb_filter_type(d1), "ZSTD")
  expect_equal(tiledb::tiledb_filter_get_option(d1, "COMPRESSION_LEVEL"), 3)
  sdf$close()
})

test_that("Metadata", {
  uri <- file.path(withr::local_tempdir(), "sdf-metadata")
  asch <- create_arrow_schema()
  sdf <- SOMADataFrameCreate(uri, asch)

  md <- list(baz = "qux", foo = "bar")
  sdf$set_metadata(md)

  # Read all metadata while the sdf is still open for write
  expect_equivalent(sdf$get_metadata("foo"), "bar")
  expect_equivalent(sdf$get_metadata("baz"), "qux")

  readmd <- sdf$get_metadata()
  expect_equivalent(readmd[["baz"]], "qux")
  expect_equivalent(readmd[["foo"]], "bar")
  sdf$close()

  # Read all metadata while the sdf is open for read
  sdf <- SOMADataFrameOpen(uri)
  readmd <- sdf$get_metadata()
  expect_equivalent(readmd[["baz"]], "qux")
  expect_equivalent(readmd[["foo"]], "bar")
  sdf$close()
})
