test_that("SOMAIndexedDataFrame creation", {
  uri <- withr::local_tempdir("soma-indexed-dataframe4")
  asch <- arrow::schema(
    foo = arrow::int32(),
    bar = arrow::float64(),
    baz = arrow::string()
  )

  sidf <- SOMAIndexedDataFrame$new(uri)
  expect_error(
    sidf$create(asch),
    "argument \"index_column_names\" is missing, with no default"
  )
  expect_error(
    sidf$create(asch, index_column_names = "qux"),
    "All 'index_column_names' must be defined in the 'schema'"
  )

  sidf$create(asch, index_column_names = "foo")
  expect_true(sidf$exists())
  expect_true(dir.exists(uri))

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

  tbl0 <- arrow::arrow_table(
    foo = 1L:10L, bar = 1.1:10.1, baz = letters[1:10]
  )
  sidf$write(tbl0)

  # read back the data (ignore attributes)
  expect_equivalent(
    tiledb::tiledb_array(sidf$uri, return_as = "asis")[],
    as.list(tbl0),
    ignore_attr = TRUE
  )

  # Read result should recreate the original Table without the soma_rowid
  tbl1 <- sidf$read()
  expect_true(tbl1$Equals(tbl0))

  # Slicing by foo
  tbl1 <- sidf$read(ids = 1L:2L)
  expect_true(tbl1$Equals(tbl1$Slice(offset = 0, length = 2)))

  # Subselecting columns
  expect_error(
    sidf$read(column_names = "foo"),
    "'column_names' must only contain non-index columns"
  )

  tbl1 <- sidf$read(column_names = "bar")
  expect_true(tbl1$Equals(tbl0$SelectColumns(c(0L, 1L))))

  # Attribute filters
  tbl1 <- sidf$read(value_filter = "bar < 5")
  expect_true(tbl1$Equals(tbl0$Filter(tbl0$bar < 5)))
})
