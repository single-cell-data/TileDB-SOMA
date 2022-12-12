test_that("SOMADataFrame creation", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- arrow::schema(
    arrow::field("foo", arrow::int32(), nullable = FALSE),
    arrow::field("bar", arrow::float64(), nullable = FALSE),
    arrow::field("baz", arrow::large_utf8(), nullable = FALSE)
  )
  sidf <- SOMADataFrame$new(uri)
  expect_error(
    sidf$create(asch),
    "argument \"index_column_names\" is missing, with no default"
  )
  expect_error(
    sidf$create(asch, index_column_names = "qux"),
    "The following field does not exist: qux"
  )

  sidf$create(asch, index_column_names = "foo")
  expect_true(sidf$exists())
  expect_true(dir.exists(uri))
  expect_match(sidf$soma_type, "SOMADataFrame")

  # check for missing columns
  expect_error(
    sidf$write(arrow::arrow_table(foo = 1L:10L)),
    "All schema fields must be present in 'values'"
  )
  # check for extra columns
  expect_error(
    sidf$write(arrow::arrow_table(qux = 1L:10L)),
    "All columns in 'values' must be defined in the schema"
  )

  tbl0 <- arrow::arrow_table(foo = 1L:10L,
                             bar = 1.1:10.1,
                             baz = letters[1:10],
                             schema=asch)

  sidf$write(tbl0)

  # read back the data (ignore attributes)
  expect_equivalent(
    tiledb::tiledb_array(sidf$uri, return_as = "asis")[],
    as.list(tbl0),
    ignore_attr = TRUE
  )

  # Read result should recreate the original Table
  tbl1 <- sidf$read()
  expect_true(tbl1$Equals(tbl0))

  # Slicing by foo
  tbl1 <- sidf$read(ids = list(foo=1L:2L))
  expect_true(tbl1$Equals(tbl1$Slice(offset = 0, length = 2)))

  # Subselecting columns
  expect_error(
    sidf$read(column_names = "foobar"),
    "'column_names' must only contain valid dimension or attribute columns"
  )

  tbl1 <- sidf$read(column_names = "bar")
  expect_true(tbl1$Equals(tbl0$SelectColumns(1L)))

  # Attribute filters
  tbl1 <- sidf$read(value_filter = "bar < 5")
  expect_true(tbl1$Equals(tbl0$Filter(tbl0$bar < 5)))
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

    ids <- bit64::as.integer64(seq(100, 109))
    sdf <- SOMADataFrame$new(uri)
    z <- sdf$read(ids = list(soma_joinid=ids))
    expect_equal(z$num_rows, 10L)
})
