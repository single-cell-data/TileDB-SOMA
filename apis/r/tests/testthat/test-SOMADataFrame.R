test_that("SOMADataFrame creation", {
  uri <- withr::local_tempdir("soma-dataframe")
  asch <- arrow::schema(
    foo = arrow::int32(),
    bar = arrow::float64(),
    baz = arrow::string()
  )

  sdf <- SOMADataFrame$new(uri)
  sdf$create(asch)

  expect_true(sdf$exists())
  expect_true(dir.exists(uri))

  tbl0 <- arrow::arrow_table(
    foo = 1L:10L, bar = 1.1:10.1, baz = letters[1:10]
  )
  expect_error(sdf$write(tbl0), "must contain a 'soma_rowid' column name")

  # add a soma_rowid column and try again
  tbl1 <- cbind(soma_rowid = 0L:9L, tbl0)
  sdf$write(tbl1)

  # read back the data (ignore attributes)
  expect_equivalent(
    tiledb::tiledb_array(sdf$uri, return_as = "asis")[],
    as.list(tbl1),
    ignore_attr = TRUE
  )

  # Read result should recreate the original Table without the soma_rowid
  tbl2 <- sdf$read()
  expect_true(tbl2$Equals(tbl0))

  # Slicing by soma_rowid
  tbl2 <- sdf$read(ids = 0:2)
  expect_true(tbl2$Equals(tbl0$Slice(offset = 0, length = 3)))

  # Subselecting columns
  tbl2 <- sdf$read(column_names = "foo")
  expect_true(tbl2$Equals(tbl0$SelectColumns(0L)))

  # Attribute filters
  tbl2 <- sdf$read(value_filter = "foo < 5")
  expect_true(tbl2$Equals(tbl0$Filter(tbl0$foo < 5)))

  # Result ordering
  expect_error(sdf$read(result_order = "foo"))
})
