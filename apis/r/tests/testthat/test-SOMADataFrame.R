test_that("SOMADataFrame creation", {
  uri <- withr::local_tempdir("soma-dataframe6")
  asch <- arrow::schema(
    foo = arrow::int32(),
    bar = arrow::float64(),
    baz = arrow::string()
  )

  sdf <- SOMADataFrame$new(uri)
  sdf$create(asch)

  expect_true(sdf$exists())
  expect_true(dir.exists(uri))

  rb0 <- arrow::record_batch(
    data.frame(foo = 1L:10L, bar = 1.1:10.1, baz = letters[1:10]),
  )
  expect_error(sdf$write(rb0), "must contain a 'soma_rowid' column name")

  # add a soma_rowid column and try again
  rb1 <- cbind(soma_rowid = 0L:9L, rb0)
  sdf$write(rb1)

  # read back the data (ignore attributes)
  expect_equivalent(
    tiledb::tiledb_array(sdf$uri, return_as = "asis")[],
    as.list(rb1),
    ignore_attr = TRUE
  )

  # Read result should recreate the original RecordBatch without the soma_rowid
  rb2 <- sdf$read()
  expect_true(rb2$Equals(rb0))

  # Slicing by soma_rowid
  rb2 <- sdf$read(ids = 0:2)
  expect_true(rb2$Equals(rb0$Slice(offset = 0, length = 3)))

  # Subselecting columns
  rb2 <- sdf$read(column_names = "foo")
  expect_true(rb2$Equals(rb0$SelectColumns(0L)))

  # Attribute filters
  rb2 <- sdf$read(value_filter = "foo < 5")
  expect_true(rb2$Equals(rb0$Filter(rb0$foo < 5)))

  # Result ordering
  expect_error(sdf$read(result_order = "foo"))
})
