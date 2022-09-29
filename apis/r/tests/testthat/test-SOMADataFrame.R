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

  rb <- arrow::record_batch(
    data.frame(foo = 1L:10L, bar = 1.1:10.1, baz = letters[1:10]),
  )
  expect_error(sdf$write(rb), "must contain a 'soma_rowid' column name")

  # add a soma_rowid column and try again
  rb <- cbind(soma_rowid = 0L:9L, rb)
  sdf$write(rb)

  # read back the data (ignore attributes)
  expect_equivalent(
    tiledb::tiledb_array(sdf$uri, return_as = "asis")[],
    as.list(rb),
    ignore_attr = TRUE
  )
})
