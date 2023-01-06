test_that("Basic mechanics", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- create_arrow_schema()

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
                             soma_joinid = 1L:10L,
                             bar = 1.1:10.1,
                             baz = letters[1:10],
                             schema = asch)

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
  tbl1 <- sidf$read(ids = list(foo = 1L:2L))
  expect_true(tbl1$Equals(tbl0$Slice(offset = 0, length = 2)))

  # Subselecting columns
  expect_error(
    sidf$read(column_names = "foobar"),
    "'column_names' must only contain valid dimension or attribute columns"
  )

  tbl1 <- sidf$read(column_names = "bar")
  expect_true(tbl1$Equals(tbl0$SelectColumns(2L)))

  # Attribute filters
  tbl1 <- sidf$read(value_filter = "bar < 5")
  expect_true(tbl1$Equals(tbl0$Filter(tbl0$bar < 5)))
})

test_that("creation with all supported dimension data types", {
  tbl0 <- arrow::arrow_table(
    double = 1.1:10.1,
    int = 1L:10L,
    int64 = bit64::as.integer64(1L:10L),
    string = letters[1:10]
  )

  for (dtype in tbl0$ColumnNames()) {
    uri <- withr::local_tempdir(paste0("soma-dataframe-", dtype))
    sidf <- SOMADataFrame$new(uri)
    expect_silent(
      sidf$create(tbl0$schema, index_column_names = dtype)
    )
    expect_true(sidf$exists())
  }
})

test_that("int64 values are stored correctly", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- arrow::schema(
    arrow::field("foo", arrow::int32(), nullable = FALSE),
    arrow::field("soma_joinid", arrow::int64(), nullable = FALSE),
  )

  sidf <- SOMADataFrame$new(uri)
  sidf$create(asch, index_column_names = "foo")
  tbl0 <- arrow::arrow_table(foo = 1L:10L, soma_joinid = 1L:10L, schema = asch)

  orig_downcast_value <- getOption("arrow.int64_downcast")

  sidf$write(tbl0)
  tbl1 <- sidf$read()
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

    ids <- bit64::as.integer64(seq(100, 109))
    sdf <- SOMADataFrame$new(uri)
    z <- sdf$read(ids = list(soma_joinid=ids))
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

  sidf <- SOMADataFrame$new(uri)
  expect_error(
    sidf$create(asch, index_column_names = "foo"),
    "Column names must not start with reserved prefix 'soma_'"
  )
})

test_that("soma_joinid is added on creation", {
  uri <- withr::local_tempdir("soma-indexed-dataframe")
  asch <- create_arrow_schema()
  asch <- asch$RemoveField(match("soma_joinid", asch$names) - 1)

  sidf <- SOMADataFrame$new(uri)
  sidf$create(asch, index_column_names = "foo")

  expect_true("soma_joinid" %in% sidf$attrnames())
  expect_equal(tiledb::datatype(sidf$attributes()$soma_joinid), "INT64")
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

  sidf <- SOMADataFrame$new(uri)
  expect_error(
    sidf$create(asch, index_column_names = "foo"),
    "soma_joinid field must be of type Arrow int64"
  )
})
