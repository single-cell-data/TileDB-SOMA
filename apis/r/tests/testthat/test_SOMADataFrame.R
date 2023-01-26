
test_that("Basic mechanics", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- create_arrow_schema()

  sdf <- SOMADataFrame$new(uri)
  expect_error(
    sdf$create(asch),
    "argument \"index_column_names\" is missing, with no default"
  )
  expect_error(
    sdf$create(asch, index_column_names = "qux"),
    "The following field does not exist: qux"
  )

  sdf$create(asch, index_column_names = "foo")
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

  tbl0 <- arrow::arrow_table(foo = 1L:36L,
                             soma_joinid = 1L:36L,
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

  # Read result should recreate the original Table
  tbl1 <- sdf$read()
  expect_true(tbl1$Equals(tbl0))

  # Slicing by foo
  tbl1 <- sdf$read(coords = list(foo = 1L:2L))
  expect_true(tbl1$Equals(tbl0$Slice(offset = 0, length = 2)))

  # Subselecting columns
  expect_error(
    sdf$read(column_names = "foobar"),
    "'column_names' must only contain valid dimension or attribute columns"
  )

  tbl1 <- sdf$read(column_names = "bar")
  expect_true(tbl1$Equals(tbl0$SelectColumns(2L)))

  # Attribute filters
  tbl1 <- sdf$read(value_filter = "bar < 5")
  expect_true(tbl1$Equals(tbl0$Filter(tbl0$bar < 5)))
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
    sdf <- SOMADataFrame$new(uri)
    expect_silent(
      sdf$create(tbl0$schema, index_column_names = dtype)
    )
    expect_true(sdf$exists())
  }
})

test_that("int64 values are stored correctly", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- arrow::schema(
    arrow::field("foo", arrow::int32(), nullable = FALSE),
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
  )

  sdf <- SOMADataFrame$new(uri)
  sdf$create(asch, index_column_names = "foo")
  tbl0 <- arrow::arrow_table(foo = 1L:10L, soma_joinid = 1L:10L, schema = asch)

  orig_downcast_value <- getOption("arrow.int64_downcast")

  sdf$write(tbl0)
  tbl1 <- sdf$read()
  expect_true(tbl1$Equals(tbl0))

  # verify int64_downcast option was restored
  expect_equal(getOption("arrow.int64_downcast"), orig_downcast_value)
})

test_that("SOMADataFrame read", {
    tdir <- tempfile()
    tgzfile <- system.file("raw-data", "soco-pbmc3k_processed-obs.tar.gz", package="tiledbsoma")
    untar(tarfile = tgzfile, exdir = tdir)

    uri <- file.path(tdir, "obs")

    sdf <- SOMADataFrame$new(uri)
    z <- sdf$read()
    expect_equal(z$num_rows, 2638L)
    expect_equal(z$num_columns, 6L)

    columns <- c("n_counts", "n_genes", "louvain")
    sdf <- SOMADataFrame$new(uri)
    z <- sdf$read(column_names=columns)
    expect_equal(z$num_columns, 3L)
    expect_equal(z$ColumnNames(), columns)

    columns <- c("n_counts", "does_not_exist")
    sdf <- SOMADataFrame$new(uri)
    expect_error(sdf$read(column_names=columns))

    coords <- bit64::as.integer64(seq(100, 109))
    sdf <- SOMADataFrame$new(uri)
    z <- sdf$read(coords = list(soma_joinid=coords))
    expect_equal(z$num_rows, 10L)
})

test_that("soma_ prefix is reserved", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- create_arrow_schema()

  # Add a soma_joinid column with the wrong type
  asch <- asch$AddField(
    i = 1,
    field = arrow::field("soma_foo", arrow::int32(), nullable = FALSE)
  )

  sdf <- SOMADataFrame$new(uri)
  expect_error(
    sdf$create(asch, index_column_names = "foo"),
    "Column names must not start with reserved prefix 'soma_'"
  )
})

test_that("soma_joinid is added on creation", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- create_arrow_schema()
  asch <- asch$RemoveField(match("soma_joinid", asch$names) - 1)

  sdf <- SOMADataFrame$new(uri)
  sdf$create(asch, index_column_names = "foo")

  expect_true("soma_joinid" %in% sdf$attrnames())
  expect_equal(tiledb::datatype(sdf$attributes()$soma_joinid), "INT64")
})

test_that("soma_joinid validations", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- create_arrow_schema()

  # Add a soma_joinid column with the wrong type
  asch <- asch$RemoveField(match("soma_joinid", asch$names) - 1)
  asch <- asch$AddField(
    i = 1,
    field = arrow::field("soma_joinid", arrow::int32(), nullable = FALSE)
  )

  sdf <- SOMADataFrame$new(uri)
  expect_error(
    sdf$create(asch, index_column_names = "foo"),
    "soma_joinid field must be of type Arrow int64"
  )
})
